#pragma once

#include "core/lp_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::simplex {

// Phase 1.1 throwaway: a dense-tableau, bounded-variable, two-phase primal
// simplex. Its only job is to prove the MPS reader, the independent
// checker and the bench harness work end to end (Checkpoint 1.1:
// `./solver netlib/afiro.mps` returns the published optimum). It rebuilds
// and inverts the dense basis matrix from scratch every iteration
// (O(rows^3) per pivot) and is expected to be too slow for anything beyond
// small Netlib instances. It is replaced by the sparse revised simplex with
// Forrest-Tomlin updates in Phase 1.2 / 2.1 and should not be extended.
core::Solution SolveDense(const core::LpProblem& problem, int max_iterations = -1,
                           const core::TolerancePolicy& tol = core::DefaultTolerances());

}  // namespace inferno::simplex
