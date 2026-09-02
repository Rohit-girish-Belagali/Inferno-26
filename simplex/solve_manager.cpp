#include "simplex/solve_manager.hpp"

#include "checker/checker.hpp"
#include "simplex/dual_simplex.hpp"
#include "simplex/revised_simplex.hpp"

namespace inferno::simplex {

ManagedSolution SolveManaged(const core::LpProblem& problem, const core::TolerancePolicy& tol) {
  core::Solution dual = SolveDual(problem, -1, tol);

  // Only ever hand back dual's answer for a checker-verified optimal.
  // dual_simplex.cpp verifies both primal and dual feasibility itself
  // before returning kOptimal (see its MaxBoundViolation/
  // MaxDualInfeasibility checks), but re-checking independently here costs
  // nothing and matches how every solve path in this project is treated —
  // trust nothing without the checker, including your own internal
  // checks. dual's kInfeasible conclusion specifically is NOT re-verified
  // this way (that path doesn't have an equivalent independent check yet
  // — a known gap, not silently relied on), so any non-optimal dual
  // status, including kInfeasible, falls through to the primal solve
  // rather than being trusted directly.
  if (dual.status == core::SolveStatus::kOptimal) {
    checker::CheckResult check = checker::VerifySolution(problem, dual);
    if (check.passed) return {dual, "dual"};
  }

  core::Solution primal = SolveRevised(problem, -1, tol);
  return {primal, "revised"};
}

}  // namespace inferno::simplex
