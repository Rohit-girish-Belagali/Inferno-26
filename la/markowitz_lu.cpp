#include "la/markowitz_lu.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace inferno::la {

bool FactorizeMarkowitz(const core::CscMatrix& b, const MarkowitzOptions& opts,
                         double singularity_tol, LuFactors& out) {
  int m = b.rows;
  if (b.cols != m) return false;
  if (m == 0) {
    out = LuFactors{};
    out.m = 0;
    return true;
  }

  // Active working matrix, kept in sync both by original column and by
  // original row so elimination can look up "which columns does the pivot
  // row touch" (needed to know what to update) as cheaply as "which rows
  // does the pivot column touch" (needed to find multipliers).
  std::vector<std::unordered_map<int, double>> active_col(m);
  std::vector<std::unordered_map<int, double>> active_row(m);
  for (int c = 0; c < m; ++c) {
    for (int p = b.col_ptr[c]; p < b.col_ptr[c + 1]; ++p) {
      int r = b.row_idx[p];
      double v = b.values[p];
      if (v == 0.0) continue;
      active_col[c][r] = v;
      active_row[r][c] = v;
    }
  }
  std::vector<char> col_active(m, 1);

  std::vector<int> row_perm(m, -1), col_perm(m, -1);
  std::vector<int> row_step(m, -1), col_step(m, -1);  // original index -> k, filled in as pivoted

  // Entries recorded during elimination, still keyed by ORIGINAL row/column
  // index (their k-index isn't known until that row/column is itself
  // pivoted, possibly at a later step) — remapped to k-space at the end.
  std::vector<std::vector<std::pair<int, double>>> l_col_orig(m);
  std::vector<std::vector<std::pair<int, double>>> u_row_orig(m);
  std::vector<double> diag_orig(m, 0.0);  // keyed by k directly (pivot found this step)

  for (int k = 0; k < m; ++k) {
    // --- Pivot search: minimum Markowitz count among threshold-stable
    // candidates, scanning every still-active column. ---
    int best_c = -1, best_r = -1;
    long long best_count = -1;

    for (int c = 0; c < m; ++c) {
      if (!col_active[c]) continue;
      const auto& colmap = active_col[c];
      if (colmap.empty()) return false;  // structurally singular

      double max_abs = 0.0;
      for (const auto& [r, v] : colmap) max_abs = std::max(max_abs, std::abs(v));
      if (max_abs <= 0.0) return false;

      for (const auto& [r, v] : colmap) {
        if (std::abs(v) < opts.pivot_threshold * max_abs) continue;
        long long count = static_cast<long long>(colmap.size() - 1) *
                           static_cast<long long>(active_row[r].size() - 1);
        if (best_c == -1 || count < best_count) {
          best_count = count;
          best_c = c;
          best_r = r;
        }
      }
    }
    if (best_c == -1) return false;

    int c = best_c, r = best_r;
    double piv = active_col[c].at(r);
    if (std::abs(piv) < singularity_tol) return false;

    // Snapshot state needed for this step before mutating it.
    std::vector<std::pair<int, double>> pivot_row_entries(active_row[r].begin(), active_row[r].end());
    std::vector<std::pair<int, double>> pivot_col_entries;
    pivot_col_entries.reserve(active_col[c].size());
    for (const auto& [i, v] : active_col[c]) {
      if (i != r) pivot_col_entries.emplace_back(i, v);
    }

    // --- Eliminate: for every other active row i with a nonzero in the
    // pivot column, subtract the appropriate multiple of the pivot row. ---
    for (const auto& [i, val_ic] : pivot_col_entries) {
      double mult = val_ic / piv;
      l_col_orig[k].emplace_back(i, mult);

      for (const auto& [cp, val_rcp] : pivot_row_entries) {
        if (cp == c) continue;
        auto& col_cp = active_col[cp];
        auto it = col_cp.find(i);
        double existing = (it != col_cp.end()) ? it->second : 0.0;
        double subtrahend = mult * val_rcp;
        double updated = existing - subtrahend;
        // Drop only genuine cancellation-to-zero, relative to the operands
        // that produced it — never an absolute threshold. On an unscaled
        // matrix an absolute cutoff (this used to be a flat 1e-13) can
        // discard a value that's small but real, and after enough
        // elimination steps that silently turns a column structurally
        // empty: FactorizeMarkowitz then reports the matrix singular when
        // it isn't. Relative-to-operands only catches true cancellation.
        double scale = std::max({std::abs(existing), std::abs(subtrahend), 1.0});
        if (std::abs(updated) < 1e-14 * scale) {
          if (it != col_cp.end()) col_cp.erase(it);
          active_row[i].erase(cp);
        } else {
          col_cp[i] = updated;
          active_row[i][cp] = updated;
        }
      }
      active_col[c].erase(i);
      active_row[i].erase(c);
    }

    // Pivot row's remaining entries (as they stood before this step's
    // elimination, which never modifies the pivot row itself) become U's
    // row k; the entry at column c is the diagonal.
    diag_orig[k] = piv;
    for (const auto& [cp, val] : pivot_row_entries) {
      if (cp == c) continue;
      u_row_orig[k].emplace_back(cp, val);
    }

    // Retire row r and column c from the active matrix.
    for (const auto& [cp, val] : pivot_row_entries) active_col[cp].erase(r);
    active_row[r].clear();
    active_col[c].clear();
    col_active[c] = 0;

    row_perm[k] = r;
    col_perm[k] = c;
    row_step[r] = k;
    col_step[c] = k;
  }

  out.m = m;
  out.row_perm = std::move(row_perm);
  out.col_perm = std::move(col_perm);
  out.row_perm_inv.assign(m, -1);
  out.col_perm_inv.assign(m, -1);
  for (int k = 0; k < m; ++k) {
    out.row_perm_inv[out.row_perm[k]] = k;
    out.col_perm_inv[out.col_perm[k]] = k;
  }

  out.diag = std::move(diag_orig);
  out.l_col.assign(m, {});
  out.u_row.assign(m, {});
  for (int k = 0; k < m; ++k) {
    for (const auto& [orig_row, mult] : l_col_orig[k]) {
      out.l_col[k].emplace_back(row_step[orig_row], mult);
    }
    for (const auto& [orig_col, val] : u_row_orig[k]) {
      out.u_row[k].emplace_back(col_step[orig_col], val);
    }
  }
  out.updates_since_refactor = 0;

  BuildTransposeViews(out);
  return true;
}

}  // namespace inferno::la
