#pragma once

#include <utility>
#include <vector>

#include "la/lu_factors.hpp"

namespace inferno::la {

// FTRAN: solves B x = b for x, given the P B Q = L U factorization. `b` is
// given sparse, as (original row index, value) pairs; returns a dense x
// indexed by original column index.
//
// Uses Gilbert & Peierls 1988's technique: a DFS over each triangular
// factor's dependency graph, seeded by b's nonzero pattern, finds the
// "reach" set — the only indices that can possibly end up nonzero — so the
// solve touches O(reach) entries instead of scanning all m unconditionally.
// This is what makes hypersparse FTRAN/BTRAN calls (a single sparse column
// against a large basis) cheap; see NOTICE_ALGORITHMS.md.
std::vector<double> Ftran(const LuFactors& lu, const std::vector<std::pair<int, double>>& b_sparse);

// BTRAN: solves B^T y = c for y, same factorization. `c` is sparse
// (original column index, value) pairs; returns dense y indexed by
// original row index.
std::vector<double> Btran(const LuFactors& lu, const std::vector<std::pair<int, double>>& c_sparse);

}  // namespace inferno::la
