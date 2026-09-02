#pragma once

#include <vector>

#include "core/sparse.hpp"

namespace inferno::la {

// Row and column scale factors for a sparse matrix: computed so that
// applying them (multiplying row i by row_scale[i] and column j by
// col_scale[j]) brings nonzero magnitudes closer to 1, which is what
// threshold pivoting in the LU factorization needs to make sound
// stability comparisons — an LP whose coefficients span 1e-6 to 1e6
// otherwise produces "small" pivots that are only small because of units,
// not because they're numerically dangerous.
struct ScaleFactors {
  std::vector<double> row_scale;
  std::vector<double> col_scale;
};

// Geometric-mean scaling (a few passes of row/column factors set to the
// inverse geometric mean of the max and min nonzero magnitude in that
// row/column under the current scaling), followed by an infinity-norm
// equilibration pass (rescale so the largest magnitude in every row and
// column is exactly 1). Citation: this is the standard two-stage scaling
// described in Suhl & Suhl 1990 ("A fast LU update for LP") and used by
// most production simplex codes; see NOTICE_ALGORITHMS.md.
ScaleFactors ComputeGeometricScaling(const core::CscMatrix& a, int geometric_passes = 2);

// Returns a new matrix with entry (i,j) multiplied by row_scale[i] *
// col_scale[j].
core::CscMatrix ApplyScaling(const core::CscMatrix& a, const ScaleFactors& scale);

}  // namespace inferno::la
