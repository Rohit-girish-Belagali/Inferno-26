#pragma once

#include "core/sparse.hpp"
#include "la/lu_factors.hpp"

namespace inferno::la {

struct MarkowitzOptions {
  // Threshold partial pivoting stability factor τ ∈ (0, 1]: a candidate
  // pivot in a column is accepted only if its magnitude is at least
  // τ * (largest magnitude in that column's currently-active part). Lower
  // τ favors sparser factors (more freedom to pick a low-fill pivot);
  // higher τ favors numerical stability. Matches
  // core::TolerancePolicy::markowitz_threshold.
  double pivot_threshold = 0.1;
};

// Factorizes a square sparse matrix B (m x m, CSC) into P B Q = L U via
// Markowitz-count pivot selection (minimize (row nnz - 1) * (col nnz - 1)
// among the active submatrix) with threshold partial pivoting for
// stability. Citation: Markowitz 1957 ("The elimination form of the
// inverse..."); threshold pivoting and the elimination bookkeeping follow
// Suhl & Suhl 1990 — see NOTICE_ALGORITHMS.md.
//
// Returns false if B is structurally or numerically singular (an active
// column becomes empty, or every remaining pivot candidate falls below
// `singularity_tol` in magnitude); `out` is left in a partially-built,
// unusable state in that case.
//
// This is a correctness-first MVP: pivot search rescans every active
// column's entries from scratch each step (no incrementally-maintained
// Markowitz counts), so it is O(m * nnz_active) rather than the O(nnz)
// amortized cost a production implementation would target. Acceptable for
// Phase 1.2's correctness gate; revisit if it becomes the bottleneck once
// wired into the real simplex in Phase 2.
bool FactorizeMarkowitz(const core::CscMatrix& b, const MarkowitzOptions& opts,
                         double singularity_tol, LuFactors& out);

}  // namespace inferno::la
