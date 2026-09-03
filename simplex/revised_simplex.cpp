#include "simplex/revised_simplex.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "core/sparse.hpp"
#include "la/basis_factorization.hpp"
#include "la/markowitz_lu.hpp"
#include "la/scaling.hpp"

namespace inferno::simplex {

namespace {

using core::kInfinity;

enum class Status { kAtLower, kAtUpper, kFree, kBasic };

// BUILD_PLAN_V2.md Phase 2.1 checklist: "cost perturbation and shifting."
// The textbook version perturbs the OBJECTIVE ITSELF, then must guarantee
// the perturbation is small enough to never change which vertex is truly
// optimal — a bound this project has no principled way to certify for an
// arbitrary Netlib instance, and getting it wrong would mean silently
// reporting a suboptimal "optimal". Bland's-rule fallback above already
// gives PROVEN, zero-risk anti-cycling (guaranteed termination, no
// wrong-answer risk); what it doesn't help with is picking a WORSE
// direction than necessary among several equally-scored candidates round
// after round on a highly degenerate problem, which is exactly what
// makes some instances grind through many more pivots than they need
// before Bland's threshold even kicks in.
//
// So: perturb only the Devex SCORE used to pick among candidates that are
// ALL improving, by a fixed, tiny, deterministic per-column factor —
// never the reduced cost itself, and never the ratio test's bound
// targets. This can change WHICH of several equally-improving directions
// gets tried first, and in what order ties are broken pivot after pivot,
// but can never change whether a direction is improving, never changes a
// ratio-test outcome, and so can never turn a correct optimum into a
// wrong one — the entering/leaving selection logic's own correctness
// conditions are completely unaffected. A deterministic hash rather than
// real randomness, so a given problem always pivots the same way run to
// run (matters for reproducing/debugging a slow instance).
double PerturbationFactor(int j) {
  uint32_t h = static_cast<uint32_t>(j) * 2654435761u;
  h ^= h >> 15;
  double frac = (h % 1000000) / 1000000.0;  // in [0, 1)
  return 1.0 + 1e-9 * frac;
}

// The column of M = [A | -I] for variable `var` in the LP's original
// row-index space — shared by FTRAN (entering column), the basis-matrix
// builder, and the eta update, so it's computed once per pivot and reused.
std::vector<std::pair<int, double>> ColumnOf(const core::LpProblem& problem, int var) {
  std::vector<std::pair<int, double>> col;
  if (var < problem.num_cols) {
    for (int p = problem.a.col_ptr[var]; p < problem.a.col_ptr[var + 1]; ++p) {
      col.emplace_back(problem.a.row_idx[p], problem.a.values[p]);
    }
  } else {
    col.emplace_back(var - problem.num_cols, -1.0);
  }
  return col;
}

core::CscMatrix BuildBasisMatrix(const core::LpProblem& problem, const std::vector<int>& basis) {
  int m = problem.num_rows;
  core::CscBuilder builder(m, m);
  for (int slot = 0; slot < m; ++slot) {
    for (const auto& [row, val] : ColumnOf(problem, basis[slot])) builder.AddEntry(slot, row, val);
  }
  return std::move(builder).Build();
}

// Crash basis (BUILD_PLAN_V2.md Phase 2.1: "do not start from the slack
// basis"): seed with structural columns instead of an all-slack start, so
// phase 1 has less work to do reaching feasibility on problems where a
// structural column is obviously a better fit for a row than its slack.
//
// Greedy weighted matching over every (row, column) entry, largest
// |coefficient| first: accept a candidate if both its row (still slack)
// and its column (not yet used) are free. This is a heuristic, not a
// guarantee of anything — the row it ends up assigned to isn't
// necessarily triangular, and the resulting basis isn't necessarily
// nonsingular. That's fine: this is purely a warm start. The caller
// always attempts to factorize whatever this returns and falls back to
// the plain all-slack basis on failure, so a bad or degenerate crash
// costs iterations, never correctness — phase 1 and phase 2 still verify
// feasibility and optimality from scratch regardless of where they
// started.
//
// A row is only crashed out of its slack when that slack has at least
// one finite bound to become nonbasic at (skips the degenerate case of a
// row with no constraint at all, row_lo=-inf and row_hi=+inf, vanishingly
// rare in practice but not worth a special-cased nonbasic-free-variable
// start).
std::vector<int> ChooseCrashBasis(const core::LpProblem& problem, const std::vector<double>& lo,
                                   const std::vector<double>& hi) {
  int m = problem.num_rows;
  std::vector<int> basis(m);
  for (int i = 0; i < m; ++i) basis[i] = problem.num_cols + i;

  struct Candidate {
    int row, col;
    double abs_val;
  };
  std::vector<Candidate> candidates;
  for (int j = 0; j < problem.num_cols; ++j) {
    for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
      int row = problem.a.row_idx[p];
      int slack = problem.num_cols + row;
      if (!std::isfinite(lo[slack]) && !std::isfinite(hi[slack])) continue;
      double v = std::abs(problem.a.values[p]);
      if (v < 1e-7) continue;  // would make a poorly-conditioned pivot; not worth it
      candidates.push_back({row, j, v});
    }
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate& a, const Candidate& b) { return a.abs_val > b.abs_val; });

  std::vector<char> row_taken(m, 0), col_used(problem.num_cols, 0);
  for (const auto& c : candidates) {
    if (row_taken[c.row] || col_used[c.col]) continue;
    row_taken[c.row] = 1;
    col_used[c.col] = 1;
    basis[c.row] = c.col;
  }
  return basis;
}

struct Workspace {
  int m = 0, n = 0;
  std::vector<double> lo, hi, cost_phase2;
  std::vector<int> basis;          // slot -> variable
  std::vector<int> basis_slot_of;  // variable -> slot, -1 if nonbasic
  std::vector<Status> status;
  std::vector<double> value;  // current value of every variable
  la::BasisFactorization bf;
  std::vector<double> devex_weight;  // Devex reference weights, meaningful for nonbasic j only
  // Set by RepairSingularBasis. A repair swaps basic columns for slacks,
  // which MOVES the current point — legitimately, but phase 2 assumes it
  // starts from the feasible point phase 1 handed it, and has no machinery
  // to restore feasibility if a mid-phase-2 repair destroys it. The driver
  // below watches this flag and restarts the two-phase sequence from the
  // repaired (valid, factorizable) basis rather than reporting a failure
  // that is really just "the point moved".
  bool repaired = false;

  // Row-major (CSR) view of A, built once per solve. The Devex weight
  // update needs the pivot ROW of the tableau against every nonbasic
  // column, which is a row-oriented question; answering it column-by-
  // column against the CSC matrix meant touching every nonzero in the
  // problem on every single pivot. See the update block for the rest.
  std::vector<int> row_ptr;
  std::vector<int> row_cols;
  std::vector<double> row_vals;
  // Scatter accumulator for that update, size n, held at all-zero between
  // pivots so it never needs an O(n) clear — only the entries actually
  // touched are reset, and `touched` records exactly which those are.
  std::vector<double> alpha_row;
  std::vector<int> touched;
};

// Builds the row-major view in `ws` from the problem's CSC matrix.
void BuildRowView(Workspace& ws, const core::LpProblem& problem) {
  int nnz = problem.a.nnz();
  ws.row_ptr.assign(ws.m + 1, 0);
  for (int p = 0; p < nnz; ++p) ++ws.row_ptr[problem.a.row_idx[p] + 1];
  for (int i = 0; i < ws.m; ++i) ws.row_ptr[i + 1] += ws.row_ptr[i];
  ws.row_cols.resize(nnz);
  ws.row_vals.resize(nnz);
  std::vector<int> pos(ws.row_ptr.begin(), ws.row_ptr.end() - 1);
  for (int j = 0; j < problem.num_cols; ++j) {
    for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
      int dst = pos[problem.a.row_idx[p]]++;
      ws.row_cols[dst] = j;
      ws.row_vals[dst] = problem.a.values[p];
    }
  }
  ws.alpha_row.assign(ws.n, 0.0);
  ws.touched.clear();
}

// Recomputes every basic variable's value from scratch via one FTRAN call,
// replacing whatever the incremental per-pivot updates had accumulated.
// The incremental update (ws.value[basis[slot]] += ... each pivot) is
// cheap but its floating-point error compounds additively over hundreds
// of pivots; dense_simplex.cpp avoids this entirely by recomputing xB from
// scratch every single iteration (part of why it's slow). This is the
// cheap middle ground: ground xB back to an accurate value whenever a
// refactorization already pays for an FTRAN-scale operation anyway,
// rather than trusting incremental drift indefinitely.
void RecomputeBasicValues(Workspace& ws, const core::LpProblem& problem) {
  std::vector<double> mx(ws.m, 0.0);
  // Walks the CSC arrays directly rather than through ColumnOf, which
  // returns a fresh std::vector by value — one heap allocation per
  // nonbasic column, every refactorization. Same traversal, no allocation.
  for (int j = 0; j < ws.n; ++j) {
    if (ws.basis_slot_of[j] != -1) continue;
    double xj = ws.value[j];
    if (xj == 0.0) continue;
    if (j < problem.num_cols) {
      for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
        mx[problem.a.row_idx[p]] += problem.a.values[p] * xj;
      }
    } else {
      mx[j - problem.num_cols] -= xj;  // slack column is -e_row
    }
  }
  std::vector<std::pair<int, double>> rhs_sparse;
  for (int row = 0; row < ws.m; ++row) {
    if (mx[row] != 0.0) rhs_sparse.emplace_back(row, -mx[row]);
  }
  std::vector<double> xb = ws.bf.Ftran(rhs_sparse);
  for (int slot = 0; slot < ws.m; ++slot) ws.value[ws.basis[slot]] = xb[slot];
}

// Singularity repair. When refactorizing the current basis fails, the LU
// reports exactly which basis columns it could not pivot and which rows
// were left uncovered (la::SingularityInfo). Swapping each of those
// columns for the slack of an uncovered row yields a basis that is
// guaranteed nonsingular: the part the LU already pivoted is
// triangularizable by construction, and the unit replacement columns
// cover precisely the missing rows.
//
// This matters because a failed refactorization was, until now, the END
// of the solve — all three of this project's remaining NUMERICAL_ERROR
// instances (scsd1, scsd8, pilot87) died at exactly this one line, mid-
// solve, after thousands of otherwise-fine iterations. Repairing and
// continuing is what production simplex codes do; abandoning a solve
// because ONE basis went singular throws away all the work done so far
// for no good reason. Correctness is unaffected either way: the repaired
// basis is just a different starting point, and phase 1 / phase 2 still
// prove feasibility and optimality from scratch, with the independent
// checker having the final say regardless.
//
// A displaced basic variable becomes nonbasic at one of its own bounds
// (or at 0 if free). Returns false only if a needed slack is somehow
// already basic — in principle unreachable, since a basic slack is a
// singleton column with Markowitz count 0 and so would have been pivoted
// first, covering its row; checked anyway rather than assumed.
bool RepairSingularBasis(Workspace& ws, const core::LpProblem& problem,
                          const la::SingularityInfo& info) {
  size_t count = std::min(info.unpivoted_cols.size(), info.unpivoted_rows.size());
  if (count == 0) return false;
  for (size_t k = 0; k < count; ++k) {
    int slot = info.unpivoted_cols[k];
    int row = info.unpivoted_rows[k];
    if (slot < 0 || slot >= ws.m || row < 0 || row >= ws.m) return false;
    int slack = problem.num_cols + row;
    if (ws.basis_slot_of[slack] != -1) return false;

    int leaving = ws.basis[slot];
    ws.basis_slot_of[leaving] = -1;
    if (std::isfinite(ws.lo[leaving])) {
      ws.status[leaving] = Status::kAtLower;
      ws.value[leaving] = ws.lo[leaving];
    } else if (std::isfinite(ws.hi[leaving])) {
      ws.status[leaving] = Status::kAtUpper;
      ws.value[leaving] = ws.hi[leaving];
    } else {
      ws.status[leaving] = Status::kFree;
      ws.value[leaving] = 0.0;
    }

    ws.basis[slot] = slack;
    ws.basis_slot_of[slack] = slot;
    ws.status[slack] = Status::kBasic;
  }
  ws.repaired = true;
  return true;
}

double MaxBoundViolation(const Workspace& ws) {
  double worst = 0.0;
  for (int j = 0; j < ws.n; ++j) {
    double v = ws.value[j];
    if (std::isfinite(ws.lo[j])) worst = std::max(worst, ws.lo[j] - v);
    if (std::isfinite(ws.hi[j])) worst = std::max(worst, v - ws.hi[j]);
  }
  return worst;
}

bool RunPhase(Workspace& ws, const core::LpProblem& problem, const std::vector<double>& cost,
              bool is_phase1, int max_iterations, const core::TolerancePolicy& tol,
              bool& numerical_error, bool& hit_iteration_limit, bool& unbounded,
              int& iterations_out) {
  numerical_error = hit_iteration_limit = unbounded = false;
  iterations_out = 0;
  int degenerate_streak = 0;
  constexpr int kBlandThreshold = 50;

  for (int iter = 0; iter < max_iterations; ++iter) {
    iterations_out = iter + 1;

    std::vector<double> cost_b(ws.m);
    if (is_phase1) {
      double total_infeas = 0.0;
      for (int slot = 0; slot < ws.m; ++slot) {
        int var = ws.basis[slot];
        double v = ws.value[var];
        if (std::isfinite(ws.lo[var]) && v < ws.lo[var] - tol.feasibility) {
          cost_b[slot] = -1.0;
          total_infeas += ws.lo[var] - v;
        } else if (std::isfinite(ws.hi[var]) && v > ws.hi[var] + tol.feasibility) {
          cost_b[slot] = 1.0;
          total_infeas += v - ws.hi[var];
        } else {
          cost_b[slot] = 0.0;
        }
      }
      if (total_infeas <= tol.feasibility) return true;
    } else {
      for (int slot = 0; slot < ws.m; ++slot) cost_b[slot] = cost[ws.basis[slot]];
    }

    std::vector<std::pair<int, double>> cost_b_sparse;
    for (int slot = 0; slot < ws.m; ++slot) {
      if (cost_b[slot] != 0.0) cost_b_sparse.emplace_back(slot, cost_b[slot]);
    }
    std::vector<double> y = ws.bf.Btran(cost_b_sparse);

    bool use_bland = degenerate_streak >= kBlandThreshold;
    int entering = -1, dir = 0;
    double best_score = 0.0;

    for (int j = 0; j < ws.n; ++j) {
      if (ws.basis_slot_of[j] != -1) continue;
      double reduced = is_phase1 ? 0.0 : cost[j];
      if (j < problem.num_cols) {
        for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
          reduced -= y[problem.a.row_idx[p]] * problem.a.values[p];
        }
      } else {
        reduced += y[j - problem.num_cols];  // Mcol_j = -e_row, so -y.Mcol = +y[row]
      }

      bool improving = false;
      int this_dir = 0;
      switch (ws.status[j]) {
        case Status::kAtLower:
          if (reduced < -tol.optimality) { improving = true; this_dir = +1; }
          break;
        case Status::kAtUpper:
          if (reduced > tol.optimality) { improving = true; this_dir = -1; }
          break;
        case Status::kFree:
          if (reduced < -tol.optimality) { improving = true; this_dir = +1; }
          else if (reduced > tol.optimality) { improving = true; this_dir = -1; }
          break;
        default:
          break;
      }
      if (!improving) continue;

      if (use_bland) {
        entering = j;
        dir = this_dir;
        break;
      }
      // Devex pricing: score by reduced^2 / reference-weight rather than
      // plain |reduced| (Dantzig). The weight approximates each nonbasic
      // direction's steepest-edge norm far more cheaply than computing it
      // exactly; see the weight-update block below the ratio test, and
      // NOTICE_ALGORITHMS.md for the citation (Harris 1973). The tiny
      // per-column perturbation factor is tie-breaking only — see
      // PerturbationFactor's own comment for why it can't affect
      // correctness.
      double score = (reduced * reduced) / ws.devex_weight[j] * PerturbationFactor(j);
      if (entering == -1 || score > best_score) {
        best_score = score;
        entering = j;
        dir = this_dir;
      }
    }

    if (entering == -1) {
      // Phase 1 legitimately ends here even while infeasible (that's how
      // the caller detects genuine infeasibility — checked separately,
      // above). Phase 2 starts from a state phase 1 already verified
      // feasible, so if phase 2 reaches "no improving direction" while
      // actually infeasible, feasibility was silently lost along the way
      // (an ill-conditioned pivot sequence, same class of risk the
      // unbounded-conclusion check above guards) — that is not a
      // trustworthy "optimal", regardless of what the reduced costs say.
      if (!is_phase1) {
        RecomputeBasicValues(ws, problem);  // rule out mere incremental drift before judging
        if (MaxBoundViolation(ws) > tol.checker_residual) {
          numerical_error = true;
          return false;
        }
      }
      return true;  // optimal for this phase
    }

    std::vector<std::pair<int, double>> entering_col = ColumnOf(problem, entering);
    std::vector<double> alpha = ws.bf.Ftran(entering_col);

    double self_limit = kInfinity;
    if (std::isfinite(ws.lo[entering]) && std::isfinite(ws.hi[entering])) {
      self_limit = ws.hi[entering] - ws.lo[entering];
    }

    // Two-pass ratio test: pass 1 finds the minimum step length; pass 2,
    // among every candidate tied with that minimum (within a small
    // tolerance), picks the numerically most stable one — largest |rate| —
    // rather than whichever was found first. A naive single-pass test can
    // pick an arbitrarily tiny pivot element on a near-tie; each such
    // choice individually clears the tol.pivot cutoff but repeatedly
    // picking the worst of several tied options measurably degrades the
    // basis over many pivots, up to and including landing on a genuinely
    // singular basis (this is the specific bug that motivated adding this
    // — see git history). Lighter than a full Harris two-pass test (no
    // bound relaxation), but fixes that failure mode.
    struct RatioCandidate {
      int slot;
      double t_limit;
      double rate_abs;
      bool to_upper;
    };
    std::vector<RatioCandidate> candidates;

    for (int slot = 0; slot < ws.m; ++slot) {
      double rate = -dir * alpha[slot];
      if (std::abs(rate) < tol.pivot) continue;
      int var = ws.basis[slot];
      double v = ws.value[var];
      double lo = ws.lo[var], hi = ws.hi[var];
      bool infeasible_below = std::isfinite(lo) && v < lo - tol.feasibility;
      bool infeasible_above = std::isfinite(hi) && v > hi + tol.feasibility;

      double t_limit = kInfinity;
      bool candidate_to_upper = false;
      if (!infeasible_below && !infeasible_above) {
        if (rate > tol.pivot && std::isfinite(hi)) {
          t_limit = (hi - v) / rate;
          candidate_to_upper = true;
        } else if (rate < -tol.pivot && std::isfinite(lo)) {
          t_limit = (lo - v) / rate;
          candidate_to_upper = false;
        }
      } else if (infeasible_below) {
        if (rate > tol.pivot) { t_limit = (lo - v) / rate; candidate_to_upper = false; }
      } else {
        if (rate < -tol.pivot) { t_limit = (hi - v) / rate; candidate_to_upper = true; }
      }
      if (!std::isfinite(t_limit)) continue;
      if (t_limit < -tol.feasibility) t_limit = 0.0;
      candidates.push_back({slot, t_limit, std::abs(rate), candidate_to_upper});
    }

    double best_t = self_limit;
    for (const auto& c : candidates) best_t = std::min(best_t, c.t_limit);

    if (!std::isfinite(best_t)) {
      // A mathematically valid "unbounded" conclusion requires standing at
      // a genuinely feasible point and finding no blocking direction — not
      // an artifact of an already-corrupted state. An ill-conditioned
      // basis (a legitimate risk with Dantzig pricing and no Harris ratio
      // test yet — see file header) can produce enormous FTRAN entries
      // that fail every ratio-test branch for a variable already grossly
      // infeasible; declaring that "unbounded" would be a confident wrong
      // answer, worse than admitting the state is untrustworthy.
      if (MaxBoundViolation(ws) > tol.checker_residual) {
        numerical_error = true;
        return false;
      }
      unbounded = true;
      return false;
    }

    double tie_band = std::max(tol.feasibility, best_t * 1e-9);
    int leaving_slot = -1;
    bool leaving_to_upper = false;
    // A bound flip never touches the basis (perfectly stable), so it wins
    // any tie against an actual pivot — give it an unbeatable score.
    double chosen_stability = (std::isfinite(self_limit) && self_limit <= best_t + tie_band)
                                   ? kInfinity
                                   : -1.0;
    for (const auto& c : candidates) {
      if (c.t_limit <= best_t + tie_band && c.rate_abs > chosen_stability) {
        chosen_stability = c.rate_abs;
        leaving_slot = c.slot;
        leaving_to_upper = c.to_upper;
      }
    }

    double t = std::max(0.0, best_t);
    if (t < tol.feasibility) ++degenerate_streak; else degenerate_streak = 0;

    for (int slot = 0; slot < ws.m; ++slot) ws.value[ws.basis[slot]] += (-dir) * alpha[slot] * t;
    ws.value[entering] += dir * t;

    if (leaving_slot == -1) {
      ws.status[entering] = (dir > 0) ? Status::kAtUpper : Status::kAtLower;
      ws.value[entering] = (dir > 0) ? ws.hi[entering] : ws.lo[entering];
    } else {
      int leaving_var = ws.basis[leaving_slot];
      ws.status[leaving_var] = leaving_to_upper ? Status::kAtUpper : Status::kAtLower;
      ws.value[leaving_var] = leaving_to_upper ? ws.hi[leaving_var] : ws.lo[leaving_var];
      ws.basis_slot_of[leaving_var] = -1;
      ws.status[entering] = Status::kBasic;
      ws.basis_slot_of[entering] = leaving_slot;
      ws.basis[leaving_slot] = entering;

      // Devex weight update. Needs the pivot row of the tableau — read it
      // off with one extra BTRAN seeded at the pivot slot (distinct from
      // the pricing BTRAN above, which is seeded by cost, not by row
      // index) — then, for every nonbasic column, how much of that row it
      // contributes. leaving_var just went nonbasic and gets the new
      // reference weight directly; every other nonbasic only grows its
      // weight, never shrinks it (matches the standard Devex update).
      {
        double alpha_rq = alpha[leaving_slot];
        std::vector<double> rho = ws.bf.Btran({{leaving_slot, 1.0}});
        double w_q = ws.devex_weight[entering];

        // alpha_r (the pivot row of the tableau) = rho^T M. Computed by
        // scattering rho's NONZEROS through the row-major view rather
        // than dotting each of the n columns against a dense rho: rho is
        // one BTRAN of a unit vector and in practice is very sparse, so
        // this touches only the nonzeros of the rows rho actually hits
        // instead of every nonzero in the problem, every pivot. The
        // previous form also built a fresh std::vector per column via
        // ColumnOf — a heap allocation per nonbasic column per pivot,
        // which profiling showed dominating the entire solve. Same
        // arithmetic, same result; only the loop order changed.
        for (int r = 0; r < ws.m; ++r) {
          double rho_r = rho[r];
          if (rho_r == 0.0) continue;
          for (int p = ws.row_ptr[r]; p < ws.row_ptr[r + 1]; ++p) {
            int j = ws.row_cols[p];
            if (ws.alpha_row[j] == 0.0) ws.touched.push_back(j);
            ws.alpha_row[j] += rho_r * ws.row_vals[p];
          }
          // Slack column for row r is -e_r, so its entry is just -rho[r].
          int slack = problem.num_cols + r;
          if (ws.alpha_row[slack] == 0.0) ws.touched.push_back(slack);
          ws.alpha_row[slack] -= rho_r;
        }

        for (int j : ws.touched) {
          double alpha_rj = ws.alpha_row[j];
          ws.alpha_row[j] = 0.0;  // leave the scratch all-zero for the next pivot
          if (alpha_rj == 0.0) continue;  // also covers a duplicate `touched` entry
          if (j == entering || ws.basis_slot_of[j] != -1) continue;
          double candidate = (alpha_rj / alpha_rq) * (alpha_rj / alpha_rq) * w_q;
          if (candidate > ws.devex_weight[j]) ws.devex_weight[j] = candidate;
        }
        ws.touched.clear();
        ws.devex_weight[leaving_var] = std::max(w_q / (alpha_rq * alpha_rq), 1.0);
      }

      bool update_ok = ws.bf.Update(leaving_slot, entering_col, tol);
      if (!update_ok || ws.bf.ShouldRefactorize()) {
        core::CscMatrix current = BuildBasisMatrix(problem, ws.basis);
        la::SingularityInfo singular;
        if (!ws.bf.Factorize(current, la::MarkowitzOptions{}, tol, &singular)) {
          // Repair the basis and try once more before giving up — see
          // RepairSingularBasis. Exactly one retry: if a basis built from
          // the LU's own report of what it could pivot still fails, the
          // trouble is not a repairable singular column.
          if (!RepairSingularBasis(ws, problem, singular)) {
            numerical_error = true;
            return false;
          }
          core::CscMatrix repaired = BuildBasisMatrix(problem, ws.basis);
          if (!ws.bf.Factorize(repaired, la::MarkowitzOptions{}, tol)) {
            numerical_error = true;
            return false;
          }
        }
        RecomputeBasicValues(ws, problem);
      }
    }
  }

  hit_iteration_limit = true;
  return false;
}

// Solves `problem` as given, with no scaling applied — SolveRevised()
// below is the public entry point and applies geometric-mean scaling
// around this.
core::Solution SolveRevisedCore(const core::LpProblem& problem, int max_iterations,
                                 const core::TolerancePolicy& tol) {
  core::Solution solution;

  Workspace ws;
  ws.m = problem.num_rows;
  ws.n = problem.num_cols + problem.num_rows;
  ws.lo.resize(ws.n);
  ws.hi.resize(ws.n);
  ws.cost_phase2.assign(ws.n, 0.0);
  ws.status.resize(ws.n);
  ws.value.assign(ws.n, 0.0);
  ws.basis.resize(ws.m);
  ws.basis_slot_of.assign(ws.n, -1);
  ws.devex_weight.assign(ws.n, 1.0);
  BuildRowView(ws, problem);

  for (int j = 0; j < problem.num_cols; ++j) {
    ws.lo[j] = problem.col_lo[j];
    ws.hi[j] = problem.col_hi[j];
    ws.cost_phase2[j] = problem.obj[j];
  }
  for (int i = 0; i < problem.num_rows; ++i) {
    int slack = problem.num_cols + i;
    ws.lo[slack] = problem.row_lo[i];
    ws.hi[slack] = problem.row_hi[i];
  }

  for (int i = 0; i < ws.m; ++i) {
    int slack = problem.num_cols + i;
    ws.basis[i] = slack;
    ws.basis_slot_of[slack] = i;
    ws.status[slack] = Status::kBasic;
  }
  for (int j = 0; j < problem.num_cols; ++j) {
    if (std::isfinite(ws.lo[j])) {
      ws.status[j] = Status::kAtLower;
      ws.value[j] = ws.lo[j];
    } else if (std::isfinite(ws.hi[j])) {
      ws.status[j] = Status::kAtUpper;
      ws.value[j] = ws.hi[j];
    } else {
      ws.status[j] = Status::kFree;
      ws.value[j] = 0.0;
    }
  }

  la::MarkowitzOptions opts;

  // Try the crash basis first; fall back to the plain all-slack basis
  // (always nonsingular — it's ±I) if it fails to factorize. Purely a
  // warm start, see ChooseCrashBasis's own comment for why this can never
  // cost correctness, only iterations.
  std::vector<int> crash_basis = ChooseCrashBasis(problem, ws.lo, ws.hi);
  core::CscMatrix crash_matrix = BuildBasisMatrix(problem, crash_basis);
  if (ws.bf.Factorize(crash_matrix, opts, tol)) {
    ws.basis = crash_basis;
    ws.basis_slot_of.assign(ws.n, -1);
    for (int slot = 0; slot < ws.m; ++slot) {
      int var = ws.basis[slot];
      ws.basis_slot_of[var] = slot;
      ws.status[var] = Status::kBasic;
    }
    for (int i = 0; i < ws.m; ++i) {
      int slack = problem.num_cols + i;
      if (ws.basis_slot_of[slack] != -1) continue;  // still basic, wasn't displaced
      // Displaced by a structural column; give it a nonbasic bound —
      // ChooseCrashBasis only displaces a slack with at least one finite
      // bound, so one of these is always available.
      if (std::isfinite(ws.lo[slack])) {
        ws.status[slack] = Status::kAtLower;
        ws.value[slack] = ws.lo[slack];
      } else {
        ws.status[slack] = Status::kAtUpper;
        ws.value[slack] = ws.hi[slack];
      }
    }
  } else {
    core::CscMatrix b0 = BuildBasisMatrix(problem, ws.basis);
    if (!ws.bf.Factorize(b0, opts, tol)) {
      solution.status = core::SolveStatus::kNumericalError;
      return solution;
    }
  }
  RecomputeBasicValues(ws, problem);

  int iter_cap = max_iterations > 0 ? max_iterations : (200 * (ws.m + ws.n) + 2000);

  bool numerical_error = false, hit_limit = false, unbounded = false;
  int iters_used = 0, phase1_iters = 0;
  std::vector<double> zero_cost(ws.n, 0.0);
  bool phase1_ok = false;
  bool p2_numerical_error = false, p2_hit_limit = false, p2_unbounded = false;
  int phase2_iters = 0;

  // Restart budget for singularity repairs. A repair leaves a valid basis
  // but a moved point, so the honest response is to redo phase 1 from
  // there rather than fail — see Workspace::repaired. Capped so a
  // pathological instance that repairs every round still terminates
  // instead of looping; each restart also costs a full phase 1, so this
  // stays small deliberately.
  constexpr int kMaxRepairRestarts = 4;
  int repair_restarts = 0;

 phase1_restart:
  ws.repaired = false;
  phase1_ok = RunPhase(ws, problem, zero_cost, /*is_phase1=*/true, iter_cap, tol,
                        numerical_error, hit_limit, unbounded, phase1_iters);
  iters_used += phase1_iters;

  if ((numerical_error || !phase1_ok) && ws.repaired && repair_restarts < kMaxRepairRestarts) {
    ++repair_restarts;
    numerical_error = hit_limit = unbounded = false;
    ws.devex_weight.assign(ws.n, 1.0);
    goto phase1_restart;
  }
  if (numerical_error) { solution.status = core::SolveStatus::kNumericalError; return solution; }
  if (hit_limit) { solution.status = core::SolveStatus::kIterationLimit; return solution; }
  if (!phase1_ok) { solution.status = core::SolveStatus::kNumericalError; return solution; }

  // Ground xB before trusting it for the feasibility check and carrying it
  // into phase 2 — phase 1 may have ended between refactorizations, so
  // whatever incremental drift built up since the last one is still live.
  RecomputeBasicValues(ws, problem);

  double total_infeas = 0.0;
  for (int j = 0; j < ws.n; ++j) {
    double v = ws.value[j];
    if (std::isfinite(ws.lo[j]) && v < ws.lo[j] - tol.feasibility) total_infeas += ws.lo[j] - v;
    if (std::isfinite(ws.hi[j]) && v > ws.hi[j] + tol.feasibility) total_infeas += v - ws.hi[j];
  }
  if (total_infeas > tol.checker_residual) {
    solution.status = core::SolveStatus::kInfeasible;
    return solution;
  }

  // Fresh Devex reference framework for phase 2 — the weights phase 1
  // built up approximate steepest-edge for the infeasibility objective,
  // not the real one.
  ws.devex_weight.assign(ws.n, 1.0);

  phase2_iters = 0;
  bool phase2_ok = RunPhase(ws, problem, ws.cost_phase2, /*is_phase1=*/false, iter_cap, tol,
                             p2_numerical_error, p2_hit_limit, p2_unbounded, phase2_iters);
  iters_used += phase2_iters;
  (void)phase2_ok;

  // Same story as phase 1: if a singularity repair moved the point out
  // from under phase 2, restoring feasibility is phase 1's job, so go do
  // that and come back rather than reporting a failure.
  if ((p2_numerical_error || p2_unbounded) && ws.repaired &&
      repair_restarts < kMaxRepairRestarts) {
    ++repair_restarts;
    p2_numerical_error = p2_hit_limit = p2_unbounded = false;
    ws.devex_weight.assign(ws.n, 1.0);
    goto phase1_restart;
  }
  if (p2_unbounded) { solution.status = core::SolveStatus::kUnbounded; return solution; }
  if (p2_numerical_error) { solution.status = core::SolveStatus::kNumericalError; return solution; }
  if (p2_hit_limit) { solution.status = core::SolveStatus::kIterationLimit; return solution; }

  // Final grounding: whatever solution we report should reflect the exact
  // current basis, not accumulated per-pivot drift since the last
  // refactorization.
  RecomputeBasicValues(ws, problem);

  std::vector<double> cost_b(ws.m);
  for (int slot = 0; slot < ws.m; ++slot) cost_b[slot] = ws.cost_phase2[ws.basis[slot]];
  std::vector<std::pair<int, double>> cost_b_sparse;
  for (int slot = 0; slot < ws.m; ++slot) {
    if (cost_b[slot] != 0.0) cost_b_sparse.emplace_back(slot, cost_b[slot]);
  }
  std::vector<double> y = ws.bf.Btran(cost_b_sparse);

  solution.status = core::SolveStatus::kOptimal;
  solution.x.assign(problem.num_cols, 0.0);
  for (int j = 0; j < problem.num_cols; ++j) solution.x[j] = ws.value[j];

  solution.row_activity.assign(problem.num_rows, 0.0);
  for (int i = 0; i < problem.num_rows; ++i) {
    solution.row_activity[i] = ws.value[problem.num_cols + i];
  }

  solution.y = y;
  solution.basis = ws.basis;

  solution.reduced_cost.assign(problem.num_cols, 0.0);
  std::vector<double> aty(problem.num_cols, 0.0);
  problem.a.TransposeMultiplyAdd(y, aty);
  for (int j = 0; j < problem.num_cols; ++j) solution.reduced_cost[j] = problem.obj[j] - aty[j];

  double obj = problem.obj_offset;
  for (int j = 0; j < problem.num_cols; ++j) obj += problem.obj[j] * solution.x[j];
  solution.objective_value = obj;
  solution.iterations = iters_used;

  return solution;
}

}  // namespace

core::Solution SolveRevised(const core::LpProblem& problem, int max_iterations,
                             const core::TolerancePolicy& tol) {
  // Geometric-mean + equilibration scaling (la/scaling.hpp), applied once
  // here rather than left in the basis matrices FactorizeMarkowitz sees.
  // An unscaled problem with widely varying coefficient magnitudes can
  // produce FTRAN/BTRAN results that are enormous purely from numerical
  // amplification (not genuine problem structure) — this was the direct
  // cause of scsd1's false "unbounded" report (a basic variable's value
  // was off by ~1e9) before that was caught defensively; scaling attacks
  // the actual root cause instead of only guarding against its symptom.
  la::ScaleFactors scale = la::ComputeGeometricScaling(problem.a);

  core::LpProblem scaled = problem;
  scaled.a = la::ApplyScaling(problem.a, scale);
  for (int i = 0; i < problem.num_rows; ++i) {
    scaled.row_lo[i] = std::isfinite(problem.row_lo[i]) ? problem.row_lo[i] * scale.row_scale[i]
                                                          : problem.row_lo[i];
    scaled.row_hi[i] = std::isfinite(problem.row_hi[i]) ? problem.row_hi[i] * scale.row_scale[i]
                                                          : problem.row_hi[i];
  }
  for (int j = 0; j < problem.num_cols; ++j) {
    scaled.col_lo[j] = std::isfinite(problem.col_lo[j]) ? problem.col_lo[j] / scale.col_scale[j]
                                                          : problem.col_lo[j];
    scaled.col_hi[j] = std::isfinite(problem.col_hi[j]) ? problem.col_hi[j] / scale.col_scale[j]
                                                          : problem.col_hi[j];
    scaled.obj[j] = problem.obj[j] * scale.col_scale[j];
  }

  core::Solution scaled_solution = SolveRevisedCore(scaled, max_iterations, tol);

  core::Solution solution;
  solution.status = scaled_solution.status;
  solution.iterations = scaled_solution.iterations;
  solution.basis = scaled_solution.basis;  // variable indices, unaffected by scaling

  if (scaled_solution.status != core::SolveStatus::kOptimal) return solution;

  solution.x.assign(problem.num_cols, 0.0);
  for (int j = 0; j < problem.num_cols; ++j) solution.x[j] = scaled_solution.x[j] * scale.col_scale[j];

  solution.row_activity.assign(problem.num_rows, 0.0);
  for (int i = 0; i < problem.num_rows; ++i) {
    solution.row_activity[i] = scaled_solution.row_activity[i] / scale.row_scale[i];
  }

  solution.y.assign(problem.num_rows, 0.0);
  for (int i = 0; i < problem.num_rows; ++i) solution.y[i] = scaled_solution.y[i] * scale.row_scale[i];

  solution.reduced_cost.assign(problem.num_cols, 0.0);
  for (int j = 0; j < problem.num_cols; ++j) {
    solution.reduced_cost[j] = scaled_solution.reduced_cost[j] / scale.col_scale[j];
  }

  double obj = problem.obj_offset;
  for (int j = 0; j < problem.num_cols; ++j) obj += problem.obj[j] * solution.x[j];
  solution.objective_value = obj;

  return solution;
}

}  // namespace inferno::simplex
