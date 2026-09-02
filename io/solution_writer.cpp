#include "io/solution_writer.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace inferno::io {

std::string StatusToString(core::SolveStatus status) {
  switch (status) {
    case core::SolveStatus::kOptimal:
      return "OPTIMAL";
    case core::SolveStatus::kInfeasible:
      return "INFEASIBLE";
    case core::SolveStatus::kUnbounded:
      return "UNBOUNDED";
    case core::SolveStatus::kIterationLimit:
      return "ITERATION_LIMIT";
    case core::SolveStatus::kNumericalError:
      return "NUMERICAL_ERROR";
  }
  return "UNKNOWN";
}

void WriteSolution(const std::string& path, const core::LpProblem& problem,
                    const core::Solution& solution) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("WriteSolution: cannot open file for writing: " + path);
  }
  out << std::setprecision(17);
  out << "Status " << StatusToString(solution.status) << "\n";
  out << "Objective " << solution.objective_value << "\n";
  out << "Iterations " << solution.iterations << "\n";
  out << "Columns " << problem.num_cols << "\n";
  for (int j = 0; j < problem.num_cols; ++j) {
    double rc = j < static_cast<int>(solution.reduced_cost.size()) ? solution.reduced_cost[j] : 0.0;
    out << problem.col_names[j] << " " << solution.x[j] << " " << rc << "\n";
  }
  out << "Rows " << problem.num_rows << "\n";
  for (int i = 0; i < problem.num_rows; ++i) {
    double dual = i < static_cast<int>(solution.y.size()) ? solution.y[i] : 0.0;
    out << problem.row_names[i] << " " << solution.row_activity[i] << " " << dual << "\n";
  }
}

}  // namespace inferno::io
