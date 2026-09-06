#pragma once

#include <vector>

#include "core/lp_problem.hpp"

namespace inferno::core {

// A mixed-integer linear program: an LP plus a flag per column saying
// which variables must take integer values. Deliberately a thin wrapper
// rather than a separate type — every MILP node solves an LP relaxation,
// so sharing LpProblem verbatim means branch-and-bound reuses the whole
// simplex, presolve and checker stack without translation.
struct MipProblem {
  LpProblem lp;
  std::vector<char> is_integer;  // size lp.num_cols; non-zero => integral
};

struct MipSolution {
  SolveStatus status = SolveStatus::kNumericalError;
  std::vector<double> x;
  double objective_value = 0.0;

  // The proven bound on the optimum: no integer solution can be better
  // than this. With the tree fully explored it equals the objective; if a
  // limit stopped the search early, the difference is the remaining gap
  // and is reported rather than hidden.
  double best_bound = 0.0;
  double gap = 0.0;  // relative, |obj - bound| / (|obj| + 1e-10)

  int nodes_explored = 0;
  bool proved_optimal = false;  // tree exhausted rather than cut short
};

}  // namespace inferno::core
