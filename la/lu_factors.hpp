#pragma once

#include <utility>
#include <vector>

namespace inferno::la {

// P B Q = L U, stored in "k-space" (pivot-step order), not original
// row/column order. Every array below is indexed by k = 0..m-1.
//
// FTRAN and BTRAN each need L and U from a different traversal direction
// (see lu_solve.cpp for the derivation of which), so both a column view and
// a row view of each factor are kept: l_col/u_row are the natural output of
// factorization (Markowitz elimination produces L one column at a time and
// U one row at a time); l_row/u_col are their transposes, built once after
// factorization (and again after every Forrest-Tomlin update).
struct LuFactors {
  int m = 0;

  std::vector<std::vector<std::pair<int, double>>> l_col;  // l_col[k]: (k', mult), k' > k
  std::vector<std::vector<std::pair<int, double>>> l_row;  // l_row[k]: (k', mult), k' < k
  std::vector<std::vector<std::pair<int, double>>> u_row;  // u_row[k]: (k', val), k' > k (diag excluded)
  std::vector<std::vector<std::pair<int, double>>> u_col;  // u_col[k]: (k', val), k' < k (diag excluded)
  std::vector<double> diag;                                // diag[k] = U[k][k], the pivot at step k

  std::vector<int> row_perm;      // row_perm[k] = original row index pivoted at step k
  std::vector<int> col_perm;      // col_perm[k] = original column index pivoted at step k
  std::vector<int> row_perm_inv;  // row_perm_inv[original row] = k
  std::vector<int> col_perm_inv;  // col_perm_inv[original col] = k

  // Number of Forrest-Tomlin eta updates applied since the last full
  // factorization; used by the refactorization policy.
  int updates_since_refactor = 0;
};

// Builds l_row and u_col from l_col and u_row. Called once right after
// Markowitz factorization and again after every basis update.
void BuildTransposeViews(LuFactors& lu);

}  // namespace inferno::la
