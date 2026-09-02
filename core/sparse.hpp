#pragma once

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace inferno::core {

// Compressed sparse column matrix. This is the canonical storage format for
// the constraint matrix A: simplex, LU and presolve all operate column-major
// because FTRAN/BTRAN and pivoting are column operations.
struct CscMatrix {
  int rows = 0;
  int cols = 0;
  std::vector<int> col_ptr;     // size cols + 1
  std::vector<int> row_idx;     // size nnz
  std::vector<double> values;   // size nnz

  int nnz() const { return static_cast<int>(values.size()); }

  void Validate() const {
    if (static_cast<int>(col_ptr.size()) != cols + 1) {
      throw std::invalid_argument("CscMatrix: col_ptr size must be cols + 1");
    }
    if (row_idx.size() != values.size()) {
      throw std::invalid_argument("CscMatrix: row_idx/values size mismatch");
    }
    for (int c = 0; c < cols; ++c) {
      if (col_ptr[c] > col_ptr[c + 1]) {
        throw std::invalid_argument("CscMatrix: col_ptr must be non-decreasing");
      }
    }
    for (int idx : row_idx) {
      if (idx < 0 || idx >= rows) {
        throw std::invalid_argument("CscMatrix: row index out of range");
      }
    }
  }

  // y += A * x  (dense x, dense y), x has length cols, y has length rows.
  void MultiplyAdd(const std::vector<double>& x, std::vector<double>& y) const {
    for (int c = 0; c < cols; ++c) {
      double xc = x[c];
      if (xc == 0.0) continue;
      for (int p = col_ptr[c]; p < col_ptr[c + 1]; ++p) {
        y[row_idx[p]] += values[p] * xc;
      }
    }
  }

  // y += A^T * x  (dense x has length rows, y has length cols).
  void TransposeMultiplyAdd(const std::vector<double>& x, std::vector<double>& y) const {
    for (int c = 0; c < cols; ++c) {
      double sum = 0.0;
      for (int p = col_ptr[c]; p < col_ptr[c + 1]; ++p) {
        sum += values[p] * x[row_idx[p]];
      }
      y[c] += sum;
    }
  }
};

// Compressed sparse row matrix. Kept alongside CSC for presolve row scans
// (e.g. singleton row detection) where row-major access avoids an O(nnz)
// transpose on every pass.
struct CsrMatrix {
  int rows = 0;
  int cols = 0;
  std::vector<int> row_ptr;     // size rows + 1
  std::vector<int> col_idx;     // size nnz
  std::vector<double> values;   // size nnz

  int nnz() const { return static_cast<int>(values.size()); }
};

// Builds a CSC matrix incrementally: caller appends (row, value) pairs one
// column at a time via AddColumn, matching how the MPS reader discovers
// entries (column-ordered).
class CscBuilder {
 public:
  CscBuilder(int rows, int cols) : rows_(rows), cols_(cols) {
    columns_.resize(cols);
  }

  void AddEntry(int col, int row, double value) {
    columns_[col].emplace_back(row, value);
  }

  CscMatrix Build() && {
    CscMatrix m;
    m.rows = rows_;
    m.cols = cols_;
    m.col_ptr.assign(cols_ + 1, 0);
    for (int c = 0; c < cols_; ++c) {
      std::sort(columns_[c].begin(), columns_[c].end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
      m.col_ptr[c + 1] = m.col_ptr[c] + static_cast<int>(columns_[c].size());
    }
    m.row_idx.resize(m.col_ptr[cols_]);
    m.values.resize(m.col_ptr[cols_]);
    for (int c = 0; c < cols_; ++c) {
      int p = m.col_ptr[c];
      for (auto& [row, value] : columns_[c]) {
        m.row_idx[p] = row;
        m.values[p] = value;
        ++p;
      }
    }
    return m;
  }

 private:
  int rows_;
  int cols_;
  std::vector<std::vector<std::pair<int, double>>> columns_;
};

}  // namespace inferno::core
