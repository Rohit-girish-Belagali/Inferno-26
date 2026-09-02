#include "simplex/revised_simplex.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/sparse.hpp"
#include "la/basis_factorization.hpp"
#include "la/markowitz_lu.hpp"
#include "la/scaling.hpp"

namespace inferno::simplex {

namespace {

using core::kInfinity;

enum class Status { kAtLower, kAtUpper, kFree, kBasic };

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

struct Workspace {
  int m = 0, n = 0;
  std::vector<double> lo, hi, cost_phase2;
  std::vector<int> basis;          // slot -> variable
  std::vector<int> basis_slot_of;  // variable -> slot, -1 if nonbasic
  std::vector<Status> status;
  std::vector<double> value;  // current value of every variable
  la::BasisFactorization bf;
  std::vector<double> devex_weight;  // Devex reference weights, meaningful for nonbasic j only
};

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
  for (int j = 0; j < ws.n; ++j) {
    if (ws.basis_slot_of[j] != -1) continue;
    double xj = ws.value[j];
    if (xj == 0.0) continue;
    for (const auto& [row, val] : ColumnOf(problem, j)) mx[row] += val * xj;
  }
  std::vector<std::pair<int, double>> rhs_sparse;
  for (int row = 0; row < ws.m; ++row) {
    if (mx[row] != 0.0) rhs_sparse.emplace_back(row, -mx[row]);
  }
  std::vector<double> xb = ws.bf.Ftran(rhs_sparse);
  for (int slot = 0; slot < ws.m; ++slot) ws.value[ws.basis[slot]] = xb[slot];
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
      // NOTICE_ALGORITHMS.md for the citation (Harris 1973).
      double score = (reduced * reduced) / ws.devex_weight[j];
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
        for (int j = 0; j < ws.n; ++j) {
          if (j == entering || ws.basis_slot_of[j] != -1) continue;
          double alpha_rj = 0.0;
          for (const auto& [row, val] : ColumnOf(problem, j)) alpha_rj += rho[row] * val;
          if (alpha_rj == 0.0) continue;
          double candidate = (alpha_rj / alpha_rq) * (alpha_rj / alpha_rq) * w_q;
          if (candidate > ws.devex_weight[j]) ws.devex_weight[j] = candidate;
        }
        ws.devex_weight[leaving_var] = std::max(w_q / (alpha_rq * alpha_rq), 1.0);
      }

      bool update_ok = ws.bf.Update(leaving_slot, entering_col, tol);
      if (!update_ok || ws.bf.ShouldRefactorize()) {
        core::CscMatrix current = BuildBasisMatrix(problem, ws.basis);
        if (!ws.bf.Factorize(current, la::MarkowitzOptions{}, tol)) {
          numerical_error = true;
          return false;
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

  // Initial basis is all-slack (B0 = -I), so xB = A * x_N directly — no
  // FTRAN needed yet.
  std::vector<double> x_nonbasic(problem.num_cols);
  for (int j = 0; j < problem.num_cols; ++j) x_nonbasic[j] = ws.value[j];
  std::vector<double> ax(ws.m, 0.0);
  problem.a.MultiplyAdd(x_nonbasic, ax);
  for (int i = 0; i < ws.m; ++i) ws.value[problem.num_cols + i] = ax[i];

  la::MarkowitzOptions opts;
  core::CscMatrix b0 = BuildBasisMatrix(problem, ws.basis);
  if (!ws.bf.Factorize(b0, opts, tol)) {
    solution.status = core::SolveStatus::kNumericalError;
    return solution;
  }

  int iter_cap = max_iterations > 0 ? max_iterations : (200 * (ws.m + ws.n) + 2000);

  bool numerical_error = false, hit_limit = false, unbounded = false;
  int iters_used = 0, phase1_iters = 0;
  std::vector<double> zero_cost(ws.n, 0.0);
  bool phase1_ok = RunPhase(ws, problem, zero_cost, /*is_phase1=*/true, iter_cap, tol,
                             numerical_error, hit_limit, unbounded, phase1_iters);
  iters_used += phase1_iters;

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

  bool p2_numerical_error = false, p2_hit_limit = false, p2_unbounded = false;
  int phase2_iters = 0;
  bool phase2_ok = RunPhase(ws, problem, ws.cost_phase2, /*is_phase1=*/false, iter_cap, tol,
                             p2_numerical_error, p2_hit_limit, p2_unbounded, phase2_iters);
  iters_used += phase2_iters;
  (void)phase2_ok;

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
