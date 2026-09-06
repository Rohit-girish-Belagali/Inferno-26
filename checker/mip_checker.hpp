#pragma once

#include <string>

#include "core/mip_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::checker {

// Independent verification of a MILP answer. Same discipline as the LP
// checker and for the same reason: the component that produced the answer
// does not get to decide whether it is correct. Everything below is
// recomputed from the raw MipProblem — the matrix, the bounds, the
// objective, the integrality flags — and none of the search's own
// intermediate values are consulted or trusted.
//
// It checks more than feasibility, because a MILP can be wrong in ways an
// LP cannot:
//   - integrality: every flagged variable is actually an integer
//   - bounds and constraints: the incumbent is genuinely feasible
//   - objective: the reported value matches the point that was returned
//   - bound sanity: the reported best_bound does not EXCEED the reported
//     objective (for a minimisation a bound above the incumbent is a
//     contradiction, and it is exactly what a broken pruning rule emits)
//   - gap: the reported gap matches the objective and bound given
//   - optimality: a claim of proved optimality is only accepted when the
//     gap is genuinely closed
struct MipCheckResult {
  bool passed = false;
  double max_integrality_violation = 0.0;
  double max_bound_violation = 0.0;
  double max_row_violation = 0.0;
  double objective_mismatch = 0.0;   // |reported - recomputed|
  double gap_mismatch = 0.0;         // |reported gap - recomputed gap|
  bool bound_contradicts_objective = false;
  bool optimality_claim_unsupported = false;
  std::string message;
};

MipCheckResult VerifyMipSolution(const core::MipProblem& problem, const core::MipSolution& solution,
                                  const core::TolerancePolicy& tol = core::DefaultTolerances());

}  // namespace inferno::checker
