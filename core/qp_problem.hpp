#pragma once

#include <string>
#include <vector>

#include "core/lp_problem.hpp"
#include "core/sparse.hpp"

namespace inferno::core {

// Convex quadratic program, in the form ADMM wants it:
//
//     min  1/2 x' P x + q' x     s.t.   row_lo <= A x <= row_hi
//
// P must be symmetric positive semidefinite for the problem to be convex;
// only its UPPER TRIANGLE is stored, since storing both halves invites
// them to disagree. Column bounds are not a separate concept here — a
// bound on x_j is just a row of A with a single 1 in column j, which is
// how the QPS reader emits them and what keeps the ADMM iteration to one
// projection instead of two.
struct QpProblem {
  std::string name;

  int num_rows = 0;  // rows of A
  int num_cols = 0;  // variables

  CscMatrix p_upper;  // upper triangle of P, num_cols x num_cols
  CscMatrix a;        // num_rows x num_cols

  std::vector<double> q;       // linear objective term, size num_cols
  std::vector<double> row_lo;  // size num_rows
  std::vector<double> row_hi;  // size num_rows
  double obj_offset = 0.0;

  std::vector<std::string> col_names;
  std::vector<std::string> row_names;
};

// What a QP solve produced, plus the KKT residuals an independent check
// recomputed — the same discipline the LP path follows: the solver's own
// claim of optimality is never the thing that is reported.
struct QpSolution {
  SolveStatus status = SolveStatus::kNumericalError;
  std::vector<double> x;
  std::vector<double> y;  // dual for the row constraints
  double objective_value = 0.0;
  int iterations = 0;

  double primal_residual = 0.0;    // max violation of row_lo <= Ax <= row_hi
  double dual_residual = 0.0;      // max |Px + q + A'y|
  double complementarity = 0.0;    // max violation of the sign/activity conditions
};

}  // namespace inferno::core
