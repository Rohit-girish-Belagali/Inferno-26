#pragma once

#include "core/lp_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::simplex {

// Bounded-variable dual simplex, on the same sparse LU spine as
// simplex/revised_simplex.hpp. Where primal simplex starts feasible and
// walks toward optimality, dual simplex starts *dual*-feasible (every
// reduced cost already has the right sign) and walks toward primal
// feasibility, which is why BUILD_PLAN_V2.md calls it "not optional" —
// it's what lets a later MILP branch-and-bound node warm-start from its
// parent's basis instead of resolving from scratch.
//
// Scope limit, honestly stated rather than silently worked around: this
// only supports problems where a *trivial* dual-feasible start exists —
// the all-slack basis with every nonbasic structural variable placed at
// whichever bound matches its objective coefficient's sign (cost >= 0 ->
// lower bound, cost < 0 -> upper bound). If a column's favorable bound is
// infinite (a one-sided bound in the wrong direction) or it's free with
// nonzero cost, no such trivial start exists, and the real fix — a
// composite/Big-M dual phase 1, mirroring what primal simplex's own
// phase 1 does for primal feasibility — is not implemented here; this
// solver reports kNumericalError rather than silently producing a wrong
// answer. See NOTICE_ALGORITHMS.md.
core::Solution SolveDual(const core::LpProblem& problem, int max_iterations = -1,
                          const core::TolerancePolicy& tol = core::DefaultTolerances());

}  // namespace inferno::simplex
