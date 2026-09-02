#pragma once

#include <string>

#include "core/lp_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::simplex {

// The "Concurrent Solve Manager" from the team's own architecture diagram
// (BUILD_PLAN_V2.md Phase 2 map / the SIH wireframe): race multiple solve
// paths, report whichever reaches a trustworthy answer.
//
// What this actually is right now: sequential, not concurrent. Real
// concurrency (std::thread, racing dual simplex against primal — or,
// eventually, PDLP on GPU per the wireframe's real intent — and taking
// whichever finishes first) means shared mutable state across threads,
// which is a real source of subtle bugs on its own; introducing that
// without dedicated time to get it right seemed like a worse trade than
// being upfront that this is a fallback chain, not a race. It tries dual
// simplex first (cheap to attempt, and fast when it applies — see
// dual_simplex.hpp's scope note), and falls back to the primal revised
// simplex whenever dual reports anything other than a checker-verified
// optimal or a genuine infeasible/unbounded conclusion. `winner` is set to
// which path actually produced the returned solution ("dual" or
// "revised"), so callers (and the dashboard) can report it honestly.
struct ManagedSolution {
  core::Solution solution;
  std::string winner;  // "dual" or "revised"
};

ManagedSolution SolveManaged(const core::LpProblem& problem,
                              const core::TolerancePolicy& tol = core::DefaultTolerances());

}  // namespace inferno::simplex
