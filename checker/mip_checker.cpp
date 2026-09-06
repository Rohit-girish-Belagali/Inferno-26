#include "checker/mip_checker.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

namespace inferno::checker {

using core::kInfinity;

MipCheckResult VerifyMipSolution(const core::MipProblem& problem, const core::MipSolution& solution,
                                  const core::TolerancePolicy& tol) {
  MipCheckResult r;
  const auto& lp = problem.lp;
  const int n = lp.num_cols;

  if (solution.status != core::SolveStatus::kOptimal) {
    r.passed = true;  // nothing claimed, nothing to verify
    r.message = "no optimal solution claimed; nothing to verify";
    return r;
  }
  if (static_cast<int>(solution.x.size()) != n) {
    r.message = "solution vector length does not match the problem";
    return r;
  }

  // --- Integrality, recomputed from the flags, not from the search. ---
  for (int j = 0; j < n; ++j) {
    if (!problem.is_integer[j]) continue;
    double d = std::abs(solution.x[j] - std::round(solution.x[j]));
    r.max_integrality_violation = std::max(r.max_integrality_violation, d);
  }

  // --- Column bounds. ---
  for (int j = 0; j < n; ++j) {
    if (std::isfinite(lp.col_lo[j])) {
      r.max_bound_violation = std::max(r.max_bound_violation, lp.col_lo[j] - solution.x[j]);
    }
    if (std::isfinite(lp.col_hi[j])) {
      r.max_bound_violation = std::max(r.max_bound_violation, solution.x[j] - lp.col_hi[j]);
    }
  }

  // --- Row activity, recomputed from the raw matrix. ---
  std::vector<double> ax(lp.num_rows, 0.0);
  lp.a.MultiplyAdd(solution.x, ax);
  for (int i = 0; i < lp.num_rows; ++i) {
    if (std::isfinite(lp.row_lo[i])) {
      r.max_row_violation = std::max(r.max_row_violation, lp.row_lo[i] - ax[i]);
    }
    if (std::isfinite(lp.row_hi[i])) {
      r.max_row_violation = std::max(r.max_row_violation, ax[i] - lp.row_hi[i]);
    }
  }

  // --- Objective, recomputed. ---
  double obj = lp.obj_offset;
  for (int j = 0; j < n; ++j) obj += lp.obj[j] * solution.x[j];
  r.objective_mismatch = std::abs(obj - solution.objective_value);

  // --- Bound sanity. For a minimisation the bound is a LOWER bound, so a
  // bound above the incumbent is a contradiction and means the pruning
  // rule discarded something it should not have. This is the check that
  // catches an unsound branch-and-bound, as opposed to merely an
  // infeasible answer. ---
  double scale = std::abs(solution.objective_value) + 1.0;
  if (solution.best_bound > solution.objective_value + 1e-6 * scale) {
    r.bound_contradicts_objective = true;
  }

  // --- Gap, recomputed from the reported objective and bound. ---
  double gap = std::abs(solution.objective_value - solution.best_bound) /
               (std::abs(solution.objective_value) + 1e-10);
  r.gap_mismatch = std::abs(gap - solution.gap);

  // --- An optimality claim is only credible with a closed gap. ---
  if (solution.proved_optimal && gap > 1e-6) r.optimality_claim_unsupported = true;

  r.passed = r.max_integrality_violation <= tol.checker_residual &&
             r.max_bound_violation <= tol.checker_residual &&
             r.max_row_violation <= tol.checker_residual &&
             r.objective_mismatch <= tol.checker_residual * scale &&
             r.gap_mismatch <= 1e-6 && !r.bound_contradicts_objective &&
             !r.optimality_claim_unsupported;

  std::ostringstream m;
  m << (r.passed ? "PASS" : "FAIL") << " integrality=" << r.max_integrality_violation
    << " bounds=" << r.max_bound_violation << " rows=" << r.max_row_violation
    << " objective=" << r.objective_mismatch << " gap=" << r.gap_mismatch;
  if (r.bound_contradicts_objective) m << " BOUND_ABOVE_OBJECTIVE";
  if (r.optimality_claim_unsupported) m << " OPTIMALITY_CLAIM_UNSUPPORTED";
  r.message = m.str();
  return r;
}

}  // namespace inferno::checker
