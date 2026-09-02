#include "presolve/presolve.hpp"

#include <cmath>

#include "core/sparse.hpp"

namespace inferno::presolve {

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

  bool changed = true;
  while (changed) {
    changed = false;
    for (int j = 0; j < num_cols; ++j) {
      if (!col_active[j]) continue;
      int nnz = problem.a.col_ptr[j + 1] - problem.a.col_ptr[j];

      bool is_fixed = std::isfinite(col_lo[j]) && std::isfinite(col_hi[j]) &&
                       (col_hi[j] - col_lo[j]) <= tol.feasibility;
      if (is_fixed) {
        double v = col_lo[j];
        for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
          int row = problem.a.row_idx[p];
          double coeff = problem.a.values[p];
          if (std::isfinite(row_lo[row])) row_lo[row] -= coeff * v;
          if (std::isfinite(row_hi[row])) row_hi[row] -= coeff * v;
        }
        obj_offset += obj[j] * v;
        col_active[j] = 0;
        result.stack.push_back({ReductionType::kFixedVariable, j, v});
        changed = true;
        continue;
      }

      if (nnz == 0) {
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
        result.stack.push_back({ReductionType::kEmptyColumn, j, v});
        changed = true;
        continue;
      }
    }
  }

  for (int i = 0; i < num_rows; ++i) {
    if (row_lo[i] > row_hi[i] + tol.feasibility) {
      result.infeasible = true;
      return result;
    }
  }

  // A row whose every column got fixed/removed now has a fixed (zero)
  // activity from its remaining (nonexistent) variables — its bounds must
  // contain 0, not just be a non-crossed range. row_lo <= row_hi alone
  // misses this: e.g. a row capping x <= 3 with x subsequently fixed at 5
  // shifts to bounds (-inf, -2], which doesn't cross, but a genuinely
  // empty row can't satisfy "activity <= -2" when its activity is fixed
  // at 0. This isn't full empty-row *reduction* (out of scope — see the
  // header comment), only the feasibility check that reduction would also
  // need, so a case like this doesn't silently produce a wrong answer.
  std::vector<int> active_nnz(num_rows, 0);
  for (int j = 0; j < num_cols; ++j) {
    if (!col_active[j]) continue;
    for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
      ++active_nnz[problem.a.row_idx[p]];
    }
  }
  for (int i = 0; i < num_rows; ++i) {
    if (active_nnz[i] > 0) continue;
    bool ok = (!std::isfinite(row_lo[i]) || row_lo[i] <= tol.feasibility) &&
              (!std::isfinite(row_hi[i]) || row_hi[i] >= -tol.feasibility);
    if (!ok) {
      result.infeasible = true;
      return result;
    }
  }

  std::vector<int> new_col_index(num_cols, -1);
  int reduced_n = 0;
  for (int j = 0; j < num_cols; ++j) {
    if (col_active[j]) new_col_index[j] = reduced_n++;
  }

  core::CscBuilder builder(num_rows, reduced_n);
  for (int j = 0; j < num_cols; ++j) {
    if (!col_active[j]) continue;
    int nj = new_col_index[j];
    for (int p = problem.a.col_ptr[j]; p < problem.a.col_ptr[j + 1]; ++p) {
      builder.AddEntry(nj, problem.a.row_idx[p], problem.a.values[p]);
    }
  }

  core::LpProblem& red = result.reduced;
  red.name = problem.name;
  red.num_rows = num_rows;
  red.num_cols = reduced_n;
  red.a = std::move(builder).Build();
  red.row_lo = row_lo;
  red.row_hi = row_hi;
  red.row_names = problem.row_names;
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
                          const core::Solution& reduced_solution) {
  core::Solution sol;
  sol.status = reduced_solution.status;
  if (reduced_solution.status != core::SolveStatus::kOptimal) return sol;

  sol.x.assign(original.num_cols, 0.0);
  for (int nj = 0; nj < static_cast<int>(result.reduced_col_to_original.size()); ++nj) {
    sol.x[result.reduced_col_to_original[nj]] = reduced_solution.x[nj];
  }
  // Reduction order doesn't actually matter for these two kinds (each just
  // sets one column's value independently of the others), but undo in
  // reverse for generality — future row-level reductions likely will.
  for (auto it = result.stack.rbegin(); it != result.stack.rend(); ++it) {
    sol.x[it->col] = it->value;
  }

  sol.row_activity.assign(original.num_rows, 0.0);
  original.a.MultiplyAdd(sol.x, sol.row_activity);

  // No row was ever removed or reinterpreted by this reduction set, so the
  // reduced solve's row duals apply unchanged to the original problem.
  sol.y = reduced_solution.y;

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
