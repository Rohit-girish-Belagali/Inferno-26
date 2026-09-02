#pragma once

#include <string>

#include "core/lp_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::checker {

// The independent solution checker. Standing rule (BUILD_PLAN_V2 #2): this
// runs on every solve, in every mode, forever. It must never share code
// with the solver that produced the solution under test — it recomputes
// everything from the raw problem data so a bug shared between "solve" and
// "check" cannot hide a wrong answer.
struct CheckResult {
  bool passed = false;
  double primal_residual = 0.0;      // max bound violation, row + column
  double dual_residual = 0.0;        // max violation of A^T y + z = c
  double complementarity_gap = 0.0;  // max |z_j * distance to nearest active bound|
  std::string message;
};

// Verifies a claimed-optimal solution against the raw problem. Only
// meaningful when solution.status == SolveStatus::kOptimal; for other
// statuses this function still runs the primal feasibility check on
// whatever `x` is present but does not evaluate duals.
CheckResult VerifySolution(const core::LpProblem& problem, const core::Solution& solution,
                            const core::TolerancePolicy& tol = core::DefaultTolerances());

}  // namespace inferno::checker
