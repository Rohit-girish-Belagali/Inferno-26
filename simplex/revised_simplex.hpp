#pragma once

#include "core/lp_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::simplex {

// Phase 2.1 — the real solver: bounded-variable primal revised simplex,
// built on la/'s sparse LU (Markowitz + Gilbert-Peierls FTRAN/BTRAN + PFI
// basis updates) instead of dense_simplex.hpp's throwaway dense tableau.
// Unlike the dense version, the basis is never refactorized from scratch
// every iteration — FTRAN gives the entering column, BTRAN gives pricing
// duals, and a pivot only costs an O(nnz) eta update.
//
// Pricing: Devex (reduced^2 / reference-weight, approximating steepest
// edge cheaply), falling back to Bland's rule after a run of degenerate
// (zero-length) pivots to guarantee termination — matches BUILD_PLAN_V2.md's
// stated pricing order ("Dantzig -> Devex -> steepest edge"); Dantzig was
// the first cut and Devex has since replaced it as the default (see
// NOTICE_ALGORITHMS.md for a measured, honest comparison — Devex is a
// large win on the instances it helps, e.g. 9x fewer iterations on
// d2q06c, but is not a strict improvement across the whole Netlib set).
// Not yet implemented, tracked as follow-up: full steepest-edge pricing,
// the Harris two-pass ratio test's bound-relaxation machinery (this uses
// a lighter two-pass test — see the ratio-test comment in the .cpp), dual
// simplex, cost perturbation, crash basis.
core::Solution SolveRevised(const core::LpProblem& problem, int max_iterations = -1,
                             const core::TolerancePolicy& tol = core::DefaultTolerances());

}  // namespace inferno::simplex
