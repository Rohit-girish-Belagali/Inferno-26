#include <chrono>
#include <iostream>
#include <string>

#include "checker/checker.hpp"
#include "core/lp_problem.hpp"
#include "io/mps_reader.hpp"
#include "io/solution_writer.hpp"
#include "simplex/dense_simplex.hpp"
#include "simplex/revised_simplex.hpp"

namespace {

void PrintUsage(const char* argv0) {
  std::cerr << "usage: " << argv0
            << " <problem.mps> [--solution out.sol] [--solver revised|dense]\n"
               "  --solver revised   (default) sparse LU + Gilbert-Peierls FTRAN/BTRAN + PFI update\n"
               "  --solver dense     Phase 1.1 throwaway dense tableau, kept for comparison\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    PrintUsage(argv[0]);
    return 2;
  }

  std::string mps_path = argv[1];
  std::string solution_path;
  std::string solver_name = "revised";
  for (int i = 2; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--solution" && i + 1 < argc) {
      solution_path = argv[++i];
    } else if (arg == "--solver" && i + 1 < argc) {
      solver_name = argv[++i];
    }
  }
  if (solver_name != "revised" && solver_name != "dense") {
    PrintUsage(argv[0]);
    return 2;
  }

  inferno::core::LpProblem problem;
  try {
    problem = inferno::io::ReadMps(mps_path);
  } catch (const std::exception& e) {
    std::cerr << "error reading '" << mps_path << "': " << e.what() << "\n";
    return 1;
  }

  auto start = std::chrono::steady_clock::now();
  inferno::core::Solution solution = solver_name == "dense"
                                          ? inferno::simplex::SolveDense(problem)
                                          : inferno::simplex::SolveRevised(problem);
  auto end = std::chrono::steady_clock::now();
  double elapsed_s = std::chrono::duration<double>(end - start).count();

  std::cout << problem.name << " solver=" << solver_name
            << " status=" << inferno::io::StatusToString(solution.status)
            << " objective=" << solution.objective_value << " iterations=" << solution.iterations
            << " time=" << elapsed_s << "s\n";

  if (solution.status == inferno::core::SolveStatus::kOptimal) {
    inferno::checker::CheckResult check = inferno::checker::VerifySolution(problem, solution);
    std::cout << "checker: " << check.message << "\n";
    if (!check.passed) {
      std::cerr << "checker FAILED — refusing to report this as a valid solve\n";
      return 3;
    }
  }

  if (!solution_path.empty()) {
    inferno::io::WriteSolution(solution_path, problem, solution);
  }

  return 0;
}
