#include "la/scaling.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace inferno::la {

namespace {
constexpr double kInf = std::numeric_limits<double>::infinity();
}

ScaleFactors ComputeGeometricScaling(const core::CscMatrix& a, int geometric_passes) {
  ScaleFactors scale;
  scale.row_scale.assign(a.rows, 1.0);
  scale.col_scale.assign(a.cols, 1.0);

  for (int pass = 0; pass < geometric_passes; ++pass) {
    // Row pass: hold col_scale fixed, set row_scale[i] to the inverse
    // geometric mean of the smallest and largest scaled magnitude in row i.
    std::vector<double> row_min(a.rows, kInf);
    std::vector<double> row_max(a.rows, 0.0);
    std::vector<char> row_has(a.rows, 0);
    for (int c = 0; c < a.cols; ++c) {
      for (int p = a.col_ptr[c]; p < a.col_ptr[c + 1]; ++p) {
        int row = a.row_idx[p];
        double mag = std::abs(a.values[p]) * scale.col_scale[c];
        if (mag <= 0.0) continue;
        row_min[row] = std::min(row_min[row], mag);
        row_max[row] = std::max(row_max[row], mag);
        row_has[row] = 1;
      }
    }
    for (int i = 0; i < a.rows; ++i) {
      if (row_has[i]) scale.row_scale[i] = 1.0 / std::sqrt(row_min[i] * row_max[i]);
    }

    // Column pass: hold the just-updated row_scale fixed, set col_scale[j]
    // the same way over column j.
    for (int c = 0; c < a.cols; ++c) {
      double col_min = kInf, col_max = 0.0;
      bool has = false;
      for (int p = a.col_ptr[c]; p < a.col_ptr[c + 1]; ++p) {
        double mag = std::abs(a.values[p]) * scale.row_scale[a.row_idx[p]];
        if (mag <= 0.0) continue;
        col_min = std::min(col_min, mag);
        col_max = std::max(col_max, mag);
        has = true;
      }
      if (has) scale.col_scale[c] = 1.0 / std::sqrt(col_min * col_max);
    }
  }

  // Equilibration: one infinity-norm pass, rows then columns, so the
  // largest magnitude in every row and column of the final scaled matrix
  // is (close to) exactly 1.
  {
    std::vector<double> row_max(a.rows, 0.0);
    for (int c = 0; c < a.cols; ++c) {
      for (int p = a.col_ptr[c]; p < a.col_ptr[c + 1]; ++p) {
        int row = a.row_idx[p];
        double mag = std::abs(a.values[p]) * scale.row_scale[row] * scale.col_scale[c];
        row_max[row] = std::max(row_max[row], mag);
      }
    }
    for (int i = 0; i < a.rows; ++i) {
      if (row_max[i] > 0.0) scale.row_scale[i] /= row_max[i];
    }

    for (int c = 0; c < a.cols; ++c) {
      double col_max = 0.0;
      for (int p = a.col_ptr[c]; p < a.col_ptr[c + 1]; ++p) {
        double mag = std::abs(a.values[p]) * scale.row_scale[a.row_idx[p]] * scale.col_scale[c];
        col_max = std::max(col_max, mag);
      }
      if (col_max > 0.0) scale.col_scale[c] /= col_max;
    }
  }

  return scale;
}

core::CscMatrix ApplyScaling(const core::CscMatrix& a, const ScaleFactors& scale) {
  core::CscMatrix out = a;
  for (int c = 0; c < a.cols; ++c) {
    for (int p = a.col_ptr[c]; p < a.col_ptr[c + 1]; ++p) {
      out.values[p] = a.values[p] * scale.row_scale[a.row_idx[p]] * scale.col_scale[c];
    }
  }
  return out;
}

}  // namespace inferno::la
