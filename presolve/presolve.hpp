#pragma once

#include <vector>

#include "core/lp_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::presolve {

// Phase 2.2 — shrink the problem before solving, then exactly reverse
// every reduction to recover a solution to the ORIGINAL problem.
//
// Scope of this MVP: only two reduction kinds, chosen specifically because
// neither one removes a row, which keeps postsolve's hardest part —
// recovering correct row duals for removed rows — out of scope entirely.
// A fixed variable's own reduced cost can be anything (checker/checker.cpp
// already exempts FX variables — lo == hi — from the complementarity
// check), and an empty column's reduced cost is trivially obj[j] itself
// (it appears in no row, so nothing to subtract). Row-level reductions
// (empty/singleton rows, forcing rows, duplicate rows, dominated columns,
// coefficient tightening) are BUILD_PLAN_V2.md Phase 2.2 checklist items
// NOT done here — they need real dual recovery through postsolve, tracked
// as follow-up, not implemented.
//
//   1. Fixed variable (col_lo == col_hi): substitute the known value into
//      every row it appears in (shifting that row's bounds by
//      -coefficient * value) and into the objective (as a constant
//      offset), then remove the column.
//   2. Empty column (appears in no row): fix it at whichever bound
//      minimizes its (zero-effect-on-feasibility) objective contribution,
//      then remove it. Skipped if the favorable bound is infinite — that
//      would only be genuinely unbounded, which the solver should detect
//      and report properly rather than presolve silently sidestepping.
//
// Reductions are applied to a fixpoint (fixing one variable can make
// another variable's row bounds collapse it to fixed too), then reduced
// rows/columns whose original bounds this pass could still tighten are
// exported as the final reduced problem.
enum class ReductionType { kFixedVariable, kEmptyColumn };

struct Reduction {
  ReductionType type;
  int col = -1;    // original column index
  double value = 0.0;  // the value the column was fixed to
};

struct PresolveResult {
  core::LpProblem reduced;
  std::vector<Reduction> stack;  // in application order; postsolve undoes in reverse
  bool infeasible = false;       // true if presolve itself proved infeasibility
  std::vector<int> reduced_col_to_original;  // reduced column index -> original column index
};

PresolveResult Presolve(const core::LpProblem& problem,
                         const core::TolerancePolicy& tol = core::DefaultTolerances());

// Expands a solution computed on `result.reduced` back into `original`'s
// column space, filling in every fixed/removed column's known value and
// leaving row duals exactly as the reduced solve produced them (valid
// because no row was ever removed or modified in a way that changes its
// dual's meaning — see the scope note above).
core::Solution Postsolve(const core::LpProblem& original, const PresolveResult& result,
                          const core::Solution& reduced_solution);

}  // namespace inferno::presolve
