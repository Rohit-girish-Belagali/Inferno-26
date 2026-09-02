#include "checker/checker.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

namespace inferno::checker {

using core::kInfinity;

namespace {

double BoundViolation(double value, double lo, double hi) {
  double v = 0.0;
  if (std::isfinite(lo)) v = std::max(v, lo - value);
  if (std::isfinite(hi)) v = std::max(v, value - hi);
  return v;
}

}  // namespace

CheckResult VerifySolution(const core::LpProblem& problem, const core::Solution& solution,
                            const core::TolerancePolicy& tol) {
  CheckResult result;

  if (static_cast<int>(solution.x.size()) != problem.num_cols) {
    result.message = "solution.x size does not match problem.num_cols";
    return result;
  }

  // --- Primal residual: recompute Ax ourselves rather than trusting the
  // solver's reported row_activity, so a bug that corrupts both cannot
  // slip through. ---
  std::vector<double> ax(problem.num_rows, 0.0);
  problem.a.MultiplyAdd(solution.x, ax);

  double primal_residual = 0.0;
  for (int i = 0; i < problem.num_rows; ++i) {
    primal_residual = std::max(primal_residual,
                                BoundViolation(ax[i], problem.row_lo[i], problem.row_hi[i]));
  }
  for (int j = 0; j < problem.num_cols; ++j) {
    primal_residual = std::max(
        primal_residual, BoundViolation(solution.x[j], problem.col_lo[j], problem.col_hi[j]));
  }
  result.primal_residual = primal_residual;

  if (solution.status != core::SolveStatus::kOptimal) {
    result.passed = primal_residual <= tol.checker_residual;
    result.message = result.passed ? "primal feasible (non-optimal status)"
                                    : "primal infeasible (non-optimal status)";
    return result;
  }

  bool has_duals = static_cast<int>(solution.y.size()) == problem.num_rows &&
                    static_cast<int>(solution.reduced_cost.size()) == problem.num_cols;
  if (!has_duals) {
    result.passed = primal_residual <= tol.checker_residual;
    result.message = result.passed ? "primal feasible; no duals supplied to check"
                                    : "primal infeasible; no duals supplied to check";
    return result;
  }

  // --- Dual residual: A^T y + z = c ---
  std::vector<double> aty(problem.num_cols, 0.0);
  problem.a.TransposeMultiplyAdd(solution.y, aty);

  double dual_residual = 0.0;
  for (int j = 0; j < problem.num_cols; ++j) {
    double stationarity = aty[j] + solution.reduced_cost[j] - problem.obj[j];
    dual_residual = std::max(dual_residual, std::abs(stationarity));
  }
  result.dual_residual = dual_residual;

  // --- Complementarity: z_j must be sign-consistent with the bound x_j
  // sits at, and zero if x_j is strictly interior. ---
  double gap = 0.0;
  for (int j = 0; j < problem.num_cols; ++j) {
    double lo = problem.col_lo[j];
    double hi = problem.col_hi[j];
    double x = solution.x[j];
    double z = solution.reduced_cost[j];

    bool at_lower = std::isfinite(lo) && (x - lo) <= tol.feasibility;
    bool at_upper = std::isfinite(hi) && (hi - x) <= tol.feasibility;

    if (at_lower && z < -tol.optimality) gap = std::max(gap, -z);
    if (at_upper && z > tol.optimality) gap = std::max(gap, z);
    if (!at_lower && !at_upper) gap = std::max(gap, std::abs(z));
  }
  result.complementarity_gap = gap;

  result.passed = primal_residual <= tol.checker_residual && dual_residual <= tol.checker_residual &&
                   gap <= tol.checker_residual;

  std::ostringstream msg;
  msg << (result.passed ? "PASS" : "FAIL") << " primal=" << primal_residual
      << " dual=" << dual_residual << " complementarity=" << gap;
  result.message = msg.str();
  return result;
}

}  // namespace inferno::checker
