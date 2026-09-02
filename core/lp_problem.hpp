#pragma once

#include <limits>
#include <string>
#include <vector>

#include "core/sparse.hpp"

namespace inferno::core {

constexpr double kInfinity = std::numeric_limits<double>::infinity();

enum class RowSense { kLessEqual, kGreaterEqual, kEqual, kRange };

// Canonical in-memory LP, produced by io/ readers and consumed by every
// solve path (dense proof-of-life simplex now; revised simplex, presolve,
// PDLP, MILP and QP later all share this representation). Row form:
//   row_lo[i] <= (A x)_i <= row_hi[i]
// Column form:
//   col_lo[j] <= x_j <= col_hi[j]
// Objective: minimize sum_j obj[j] * x_j + obj_offset
struct LpProblem {
  std::string name;

  int num_rows = 0;
  int num_cols = 0;

  CscMatrix a;  // rows x cols constraint matrix

  std::vector<double> row_lo;  // size num_rows
  std::vector<double> row_hi;  // size num_rows
  std::vector<std::string> row_names;

  std::vector<double> col_lo;  // size num_cols, default 0
  std::vector<double> col_hi;  // size num_cols, default +inf
  std::vector<double> obj;     // size num_cols
  double obj_offset = 0.0;
  std::vector<std::string> col_names;

  bool IsFeasibleBounds() const {
    for (int i = 0; i < num_rows; ++i) {
      if (row_lo[i] > row_hi[i]) return false;
    }
    for (int j = 0; j < num_cols; ++j) {
      if (col_lo[j] > col_hi[j]) return false;
    }
    return true;
  }
};

// Independent verdict on a candidate solution — produced by checker/, never
// by the solver that generated the solution. Kept in core/ because io/ and
// simplex/ both need to reference SolveStatus.
enum class SolveStatus {
  kOptimal,
  kInfeasible,
  kUnbounded,
  kIterationLimit,
  kNumericalError,
};

struct Solution {
  SolveStatus status = SolveStatus::kNumericalError;
  std::vector<double> x;          // primal values, size num_cols
  std::vector<double> row_activity;  // A x, size num_rows
  std::vector<double> y;          // dual values (row duals), size num_rows
  std::vector<double> reduced_cost;  // size num_cols
  double objective_value = 0.0;
  int iterations = 0;

  // The optimal basis: basis[i] is the variable index occupying basis slot
  // i (0..num_cols-1 for a structural column, num_cols..num_cols+num_rows-1
  // for slack row i's column), size num_rows. Only set when status ==
  // kOptimal; empty otherwise. Exists so downstream consumers (la/ tests,
  // eventually the real revised simplex) can reconstruct the real basis
  // matrix a solve actually produced, instead of only ever exercising
  // synthetic ones.
  std::vector<int> basis;
};

}  // namespace inferno::core
