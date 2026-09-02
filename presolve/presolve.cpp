#include "presolve/presolve.hpp"

#include <cmath>

#include "core/sparse.hpp"

namespace inferno::presolve {

namespace {

// One active row's terms, rebuilt fresh each fixpoint round rather than
// maintained incrementally — simpler and safe, and cheap at Netlib scale
// (a handful of rounds, each O(nnz)).
struct RowEntry {
  int col;
  double coeff;
};

}  // namespace

PresolveResult Presolve(const core::LpProblem& problem, const core::TolerancePolicy& tol) {
  PresolveResult result;
  int num_cols = problem.num_cols;
  int num_rows = problem.num_rows;

  std::vector<double> col_lo = problem.col_lo;
  std::vector<double> col_hi = problem.col_hi;
  std::vector<double> row_lo = problem.row_lo;
  std::vector<double> row_hi = problem.row_hi;
  std::vector<double> obj = problem.obj;
  double obj_offset = problem.obj_offset;
  std::vector<char> col_active(num_cols, 1);
  std::vector<char> row_active(num_rows, 1);

  bool changed = true;
  while (changed) {
    changed = false;

    // --- Column pass: fixed variables and empty columns. A column counts
    // as empty if every row it touches has already been removed, not just
    // if it was empty from the start. ---
    for (int j = 0; j < num_cols; ++j) {
      if (!col_active[j]) continue;
      int active_nnz = 0;
      for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
        if (row_active[problem.a.row_idx[p]]) ++active_nnz;
      }

      bool is_fixed = std::isfinite(col_lo[j]) && std::isfinite(col_hi[j]) &&
                       (col_hi[j] - col_lo[j]) <= tol.feasibility;
      if (is_fixed) {
        double v = col_lo[j];
        for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
          int row = problem.a.row_idx[p];
          if (!row_active[row]) continue;
          double coeff = problem.a.values[p];
          if (std::isfinite(row_lo[row])) row_lo[row] -= coeff * v;
          if (std::isfinite(row_hi[row])) row_hi[row] -= coeff * v;
        }
        obj_offset += obj[j] * v;
        col_active[j] = 0;
        result.stack.push_back({ReductionType::kFixedVariable, j, -1, v});
        changed = true;
        continue;
      }

      if (active_nnz == 0) {
        double v;
        if (obj[j] > tol.optimality) {
          if (!std::isfinite(col_lo[j])) continue;  // genuinely unbounded; let the solver see it
          v = col_lo[j];
        } else if (obj[j] < -tol.optimality) {
          if (!std::isfinite(col_hi[j])) continue;
          v = col_hi[j];
        } else if (std::isfinite(col_lo[j]) && col_lo[j] <= 0.0 &&
                   (!std::isfinite(col_hi[j]) || col_hi[j] >= 0.0)) {
          v = 0.0;
        } else if (std::isfinite(col_lo[j])) {
          v = col_lo[j];
        } else if (std::isfinite(col_hi[j])) {
          v = col_hi[j];
        } else {
          v = 0.0;
        }
        obj_offset += obj[j] * v;
        col_active[j] = 0;
        result.stack.push_back({ReductionType::kEmptyColumn, j, -1, v});
        changed = true;
        continue;
      }
    }

    if (result.infeasible) return result;
    for (int j = 0; j < num_cols; ++j) {
      if (col_active[j] && col_lo[j] > col_hi[j] + tol.feasibility) {
        result.infeasible = true;
        return result;
      }
    }

    // --- Build row-major view of the currently-active matrix, once per
    // round, for the row pass below. ---
    std::vector<std::vector<RowEntry>> row_terms(num_rows);
    for (int j = 0; j < num_cols; ++j) {
      if (!col_active[j]) continue;
      for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
        int row = problem.a.row_idx[p];
        if (!row_active[row]) continue;
        row_terms[row].push_back({j, problem.a.values[p]});
      }
    }

    // --- Row pass: redundant rows, singleton rows, free column
    // singletons (equality rows only), and general bound tightening. ---
    for (int i = 0; i < num_rows; ++i) {
      if (!row_active[i]) continue;
      const auto& terms = row_terms[i];

      // A row with zero active terms has activity FORCED to exactly 0 —
      // not "possibly outside bounds", but structurally, unconditionally
      // 0 — so this is an infeasibility check, not just a redundancy one:
      // the generic implied-range check below only ever REMOVES a row, it
      // never flags one infeasible, so the empty case needs its own
      // explicit branch or a genuinely-empty, out-of-bounds row would
      // silently pass through unresolved every round instead of being
      // caught.
      if (terms.empty()) {
        bool ok = (!std::isfinite(row_lo[i]) || row_lo[i] <= tol.feasibility) &&
                  (!std::isfinite(row_hi[i]) || row_hi[i] >= -tol.feasibility);
        if (!ok) {
          result.infeasible = true;
          return result;
        }
        row_active[i] = 0;
        result.stack.push_back({ReductionType::kRedundantRow, -1, i});
        changed = true;
        continue;
      }

      // Implied activity range from every active column's CURRENT bounds.
      double implied_min = 0.0, implied_max = 0.0;
      bool min_inf = false, max_inf = false;
      for (const auto& t : terms) {
        double lo = col_lo[t.col], hi = col_hi[t.col];
        if (t.coeff > 0.0) {
          if (!std::isfinite(lo)) min_inf = true; else implied_min += t.coeff * lo;
          if (!std::isfinite(hi)) max_inf = true; else implied_max += t.coeff * hi;
        } else if (t.coeff < 0.0) {
          if (!std::isfinite(hi)) min_inf = true; else implied_min += t.coeff * hi;
          if (!std::isfinite(lo)) max_inf = true; else implied_max += t.coeff * lo;
        }
      }

      bool lo_side_ok = !std::isfinite(row_lo[i]) ||
                         (!min_inf && implied_min >= row_lo[i] - tol.feasibility);
      bool hi_side_ok = !std::isfinite(row_hi[i]) ||
                         (!max_inf && implied_max <= row_hi[i] + tol.feasibility);
      if (lo_side_ok && hi_side_ok) {
        row_active[i] = 0;
        result.stack.push_back({ReductionType::kRedundantRow, -1, i});
        changed = true;
        continue;
      }

      if (terms.size() == 1) {
        int j = terms[0].col;
        double a = terms[0].coeff;
        double implied_lo, implied_hi;
        if (a > 0.0) {
          implied_lo = std::isfinite(row_lo[i]) ? row_lo[i] / a : -core::kInfinity;
          implied_hi = std::isfinite(row_hi[i]) ? row_hi[i] / a : core::kInfinity;
        } else {
          implied_lo = std::isfinite(row_hi[i]) ? row_hi[i] / a : -core::kInfinity;
          implied_hi = std::isfinite(row_lo[i]) ? row_lo[i] / a : core::kInfinity;
        }
        double prev_lo = col_lo[j], prev_hi = col_hi[j];
        double new_lo = std::max(prev_lo, implied_lo);
        double new_hi = std::min(prev_hi, implied_hi);
        if (new_lo > new_hi + tol.feasibility) {
          result.infeasible = true;
          return result;
        }
        col_lo[j] = new_lo;
        col_hi[j] = new_hi;
        row_active[i] = 0;
        Reduction r{ReductionType::kSingletonRow, j, i};
        r.coeff = a;
        r.prev_lo = prev_lo;
        r.prev_hi = prev_hi;
        r.new_lo = new_lo;
        r.new_hi = new_hi;
        result.stack.push_back(r);
        changed = true;
        continue;
      }

      // Free column singleton: exactly one term, that column is free
      // (unbounded both sides), and the row is an equality — so the row
      // is really just x_j's definition, uniquely solvable.
      //
      // x_j = (RHS - sum_{k != j} a_ik*x_k) / a_ij, so obj_j*x_j is a
      // LINEAR function of the other columns in this row, not a constant
      // — substitute it into their objective coefficients (and into
      // obj_offset for the constant term) before dropping column j,
      // exactly the way kFixedVariable folds a constant value into
      // obj_offset. Skipping this substitution was a real bug caught by
      // the full Netlib run (capri, wrong objective) — the reduced
      // problem's objective silently stopped being equal to the
      // original's.
      if (std::isfinite(row_lo[i]) && std::isfinite(row_hi[i]) &&
          (row_hi[i] - row_lo[i]) <= tol.feasibility) {
        for (const auto& t : terms) {
          int j = t.col;
          if (std::isfinite(col_lo[j]) || std::isfinite(col_hi[j])) continue;
          // Is j a singleton column across the WHOLE (still-active) matrix,
          // not just this row?
          int col_active_nnz = 0;
          for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
            if (row_active[problem.a.row_idx[p]]) ++col_active_nnz;
          }
          if (col_active_nnz != 1) continue;

          double obj_j = obj[j];
          double a_ij = t.coeff;
          if (obj_j != 0.0) {
            obj_offset += obj_j * row_lo[i] / a_ij;
            for (const auto& u : terms) {
              if (u.col == j) continue;
              obj[u.col] -= obj_j * u.coeff / a_ij;
            }
          }

          row_active[i] = 0;
          col_active[j] = 0;
          Reduction r{ReductionType::kFreeColumnSingleton, j, i};
          r.coeff = a_ij;
          r.value = row_lo[i];  // == row_hi[i], the equality RHS
          for (const auto& u : terms) {
            if (u.col == j) continue;
            r.other_active_terms.emplace_back(u.col, u.coeff);
          }
          result.stack.push_back(r);
          changed = true;
          break;
        }
        if (changed) continue;
      }

      // General (non-removing) bound tightening was tried and reverted:
      // it's sound for FEASIBILITY (never cuts off the true optimum —
      // verified against the full Netlib set) but breaks the checker,
      // because it narrows a column's bound below what's recorded in the
      // ORIGINAL problem without leaving any record for postsolve to
      // undo. When the tightened bound (not the original one) ends up
      // binding, the checker — which validates against the ORIGINAL
      // problem's bounds — sees x_j sitting strictly inside its true
      // bounds and demands reduced_cost_j == 0, but the reduced solve's
      // reduced cost reflects the SYNTHETIC tightened bound instead, and
      // is generally nonzero. Unlike kSingletonRow, the row responsible
      // here stays ACTIVE and typically has other columns depending on
      // its dual too, so the same "just set y_i to make this one
      // column's reduced cost 0" fix isn't available without risking
      // breaking THEIR complementarity — real dual reconciliation for
      // that case is a harder problem than this pass takes on. Tracked as
      // a real, not-yet-implemented Phase 2.2 checklist item, same as
      // forcing rows below.
    }
  }

  for (int i = 0; i < num_rows; ++i) {
    if (row_active[i] && row_lo[i] > row_hi[i] + tol.feasibility) {
      result.infeasible = true;
      return result;
    }
  }

  std::vector<int> new_col_index(num_cols, -1);
  int reduced_n = 0;
  for (int j = 0; j < num_cols; ++j) {
    if (col_active[j]) new_col_index[j] = reduced_n++;
  }
  std::vector<int> new_row_index(num_rows, -1);
  int reduced_m = 0;
  for (int i = 0; i < num_rows; ++i) {
    if (row_active[i]) new_row_index[i] = reduced_m++;
  }

  core::CscBuilder builder(reduced_m, reduced_n);
  for (int j = 0; j < num_cols; ++j) {
    if (!col_active[j]) continue;
    int nj = new_col_index[j];
    for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
      int row = problem.a.row_idx[p];
      if (!row_active[row]) continue;
      builder.AddEntry(nj, new_row_index[row], problem.a.values[p]);
    }
  }

  core::LpProblem& red = result.reduced;
  red.name = problem.name;
  red.num_rows = reduced_m;
  red.num_cols = reduced_n;
  red.a = std::move(builder).Build();
  red.row_lo.resize(reduced_m);
  red.row_hi.resize(reduced_m);
  red.row_names.resize(reduced_m);
  result.reduced_row_to_original.resize(reduced_m);
  for (int i = 0; i < num_rows; ++i) {
    if (!row_active[i]) continue;
    int ni = new_row_index[i];
    red.row_lo[ni] = row_lo[i];
    red.row_hi[ni] = row_hi[i];
    red.row_names[ni] = problem.row_names[i];
    result.reduced_row_to_original[ni] = i;
  }
  red.col_lo.resize(reduced_n);
  red.col_hi.resize(reduced_n);
  red.obj.resize(reduced_n);
  red.col_names.resize(reduced_n);
  result.reduced_col_to_original.resize(reduced_n);
  for (int j = 0; j < num_cols; ++j) {
    if (!col_active[j]) continue;
    int nj = new_col_index[j];
    red.col_lo[nj] = col_lo[j];
    red.col_hi[nj] = col_hi[j];
    red.obj[nj] = obj[j];
    red.col_names[nj] = problem.col_names[j];
    result.reduced_col_to_original[nj] = j;
  }
  red.obj_offset = obj_offset;

  return result;
}

core::Solution Postsolve(const core::LpProblem& original, const PresolveResult& result,
                          const core::Solution& reduced_solution, const core::TolerancePolicy& tol) {
  core::Solution sol;
  sol.status = reduced_solution.status;
  if (reduced_solution.status != core::SolveStatus::kOptimal) return sol;

  // --- Primal fill-in: reduced-problem columns keep their solved value;
  // every fixed/empty-column reduction sets its own value directly. Free
  // column singletons are filled in below, once row activity from the
  // OTHER columns is known. ---
  sol.x.assign(original.num_cols, 0.0);
  std::vector<char> x_known(original.num_cols, 0);
  for (int nj = 0; nj < static_cast<int>(result.reduced_col_to_original.size()); ++nj) {
    int j = result.reduced_col_to_original[nj];
    sol.x[j] = reduced_solution.x[nj];
    x_known[j] = 1;
  }
  for (const auto& r : result.stack) {
    if (r.type == ReductionType::kFixedVariable || r.type == ReductionType::kEmptyColumn) {
      sol.x[r.col] = r.value;
      x_known[r.col] = 1;
    }
  }
  // Free column singletons: x_j = (row RHS - sum of other columns' terms) / a.
  //
  // Must sum over `r.other_active_terms` — the columns still active in
  // this row AT APPLICATION TIME — not every column the ORIGINAL matrix
  // lists for this row. `r.value` (the RHS) was already shifted, at
  // presolve time, to exclude whatever was fixed in this row EARLIER
  // (kFixedVariable folds a fixed column's contribution into row bounds
  // the same way it folds into obj_offset); re-including those same
  // already-fixed columns here would double-count them. This was a real
  // bug the full Netlib run caught: on greenbeb, a 200+-term equality row
  // where most columns were already fixed by the time it became a
  // free-column-singleton candidate, double-counting them threw the
  // recovered value off by thousands — primal residual ~10^4, not a
  // rounding-scale miss.
  //
  // Processed in REVERSE stack order, not forward: row i2's "other
  // columns" (per `other_active_terms`) can include a column j3 that is
  // itself a LATER-applied (so later-in-stack) free column singleton —
  // j2's row, at the time it was eliminated, still had j3 as an ordinary
  // active column; j3 only became eligible for its own substitution
  // afterward, once row i2 was gone and its own active-row count dropped
  // to one. Filling j2's value in forward order would read j3's value
  // before it's known (the other real bug the full Netlib run caught —
  // capri and vtp.base's postsolved solutions failed the checker).
  // Reverse order resolves the later-applied j3 first, exactly like the
  // dual-recovery pass below.
  for (auto it = result.stack.rbegin(); it != result.stack.rend(); ++it) {
    const Reduction& r = *it;
    if (r.type != ReductionType::kFreeColumnSingleton) continue;
    double others = 0.0;
    for (const auto& [col, coeff] : r.other_active_terms) {
      others += coeff * sol.x[col];  // known by now: either resolved above, or never removed
    }
    sol.x[r.col] = (r.value - others) / r.coeff;
    x_known[r.col] = 1;
  }

  sol.row_activity.assign(original.num_rows, 0.0);
  original.a.MultiplyAdd(sol.x, sol.row_activity);

  // --- Dual recovery: reduced-problem rows keep the reduced solve's dual.
  // Removed rows get theirs from each reduction's own rule (see
  // presolve.hpp), but NOT via a single reverse sweep — unlike primal
  // fill-in, a removed row's y can depend on ANOTHER removed row's y in
  // either application-order direction: kSingletonRow's and
  // kFreeColumnSingleton's formulas both sum over every OTHER row the
  // column touches in the ORIGINAL matrix, and that set can include rows
  // removed earlier (a free-column-singleton whose "other" column j later
  // became a singleton-row target elsewhere) as well as rows removed
  // later — a genuine dependency DAG, not a stack. So: iterate to a
  // fixpoint, resolving whatever can be resolved (all of a reduction's
  // OTHER touched rows already known) each pass, until nothing changes.
  // Terminates because the dependency graph has no cycles (each row's y
  // is a function of the FINAL x and of already-fixed original data, not
  // of anything that could reference itself).
  std::vector<double> y(original.num_rows, 0.0);
  std::vector<char> row_resolved(original.num_rows, 0);
  for (int i = 0; i < original.num_rows; ++i) row_resolved[i] = 1;  // default: active rows
  for (int ni = 0; ni < static_cast<int>(result.reduced_row_to_original.size()); ++ni) {
    y[result.reduced_row_to_original[ni]] = reduced_solution.y[ni];
  }
  for (const auto& r : result.stack) {
    if (r.type == ReductionType::kRedundantRow || r.type == ReductionType::kSingletonRow ||
        r.type == ReductionType::kFreeColumnSingleton) {
      row_resolved[r.row] = 0;
    }
  }

  auto ready = [&](const Reduction& r) {
    for (int p = original.a.col_ptr[r.col]; p < original.a.col_ptr[r.col + 1]; ++p) {
      int row = original.a.row_idx[p];
      if (row != r.row && !row_resolved[row]) return false;
    }
    return true;
  };
  auto resolve_via_zero_reduced_cost = [&](const Reduction& r) {
    double partial_rc = original.obj[r.col];
    for (int p = original.a.col_ptr[r.col]; p < original.a.col_ptr[r.col + 1]; ++p) {
      int row = original.a.row_idx[p];
      if (row == r.row) continue;
      partial_rc -= original.a.values[p] * y[row];
    }
    y[r.row] = partial_rc / r.coeff;
  };

  bool progress = true;
  while (progress) {
    progress = false;
    for (const auto& r : result.stack) {
      if (r.type != ReductionType::kRedundantRow && r.type != ReductionType::kSingletonRow &&
          r.type != ReductionType::kFreeColumnSingleton) {
        continue;
      }
      if (row_resolved[r.row]) continue;

      // Whether this row needs the zero-reduced-cost formula (and hence
      // needs `ready()`, i.e. every OTHER row its column touches already
      // resolved) is itself always immediately decidable — it only reads
      // sol.x (already fully known) and the reduction's own recorded
      // bounds, never another row's y. Gating THAT decision on ready()
      // too was a real bug: two rows that each turn out to be "not
      // responsible" (both trivially 0, no real dependency) but that
      // happen to reference each other's row via a shared column would
      // deadlock — neither ever seen as ready — silently leaving them at
      // the initial y=0 default, which is only sometimes the right
      // answer. Only the formula branch itself is gated.
      bool needs_formula;
      if (r.type == ReductionType::kRedundantRow) {
        needs_formula = false;
      } else if (r.type == ReductionType::kFreeColumnSingleton) {
        needs_formula = true;  // free variable is always basic: reduced_cost_j must be exactly 0
      } else {  // kSingletonRow
        double x_j = sol.x[r.col];
        bool at_new_lo = std::isfinite(r.new_lo) && std::abs(x_j - r.new_lo) <= tol.feasibility;
        bool at_new_hi = std::isfinite(r.new_hi) && std::abs(x_j - r.new_hi) <= tol.feasibility;
        bool lo_tightened = !std::isfinite(r.prev_lo) || r.new_lo > r.prev_lo + tol.feasibility;
        bool hi_tightened = !std::isfinite(r.prev_hi) || r.new_hi < r.prev_hi - tol.feasibility;
        needs_formula = (at_new_lo && lo_tightened) || (at_new_hi && hi_tightened);
      }

      if (!needs_formula) {
        y[r.row] = 0.0;
      } else {
        if (!ready(r)) continue;
        resolve_via_zero_reduced_cost(r);
      }
      row_resolved[r.row] = 1;
      progress = true;
    }
    // kFixedVariable / kEmptyColumn don't touch row duals.
  }
  sol.y = y;

  sol.reduced_cost.assign(original.num_cols, 0.0);
  std::vector<double> aty(original.num_cols, 0.0);
  original.a.TransposeMultiplyAdd(sol.y, aty);
  for (int j = 0; j < original.num_cols; ++j) sol.reduced_cost[j] = original.obj[j] - aty[j];

  double obj = original.obj_offset;
  for (int j = 0; j < original.num_cols; ++j) obj += original.obj[j] * sol.x[j];
  sol.objective_value = obj;
  sol.iterations = reduced_solution.iterations;

  return sol;
}

}  // namespace inferno::presolve
