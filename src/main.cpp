#include <chrono>
#include <iostream>
#include <string>

#include "checker/checker.hpp"
#include "core/lp_problem.hpp"
#include "io/mps_reader.hpp"
#include "io/solution_writer.hpp"
#include "presolve/presolve.hpp"
#include "simplex/dense_simplex.hpp"
#include "simplex/dual_simplex.hpp"
#include "simplex/revised_simplex.hpp"
#include "simplex/solve_manager.hpp"

namespace {

void PrintUsage(const char* argv0) {
  std::cerr << "usage: " << argv0
            << " <problem.mps> [--solution out.sol] [--solver revised|dense|dual|managed] "
               "[--presolve]\n"
               "  --solver revised   (default) sparse LU + Gilbert-Peierls FTRAN/BTRAN + PFI update\n"
               "  --solver dense     Phase 1.1 throwaway dense tableau, kept for comparison\n"
               "  --solver dual      bounded-variable dual simplex; only when a trivial "
               "dual-feasible start exists (see simplex/dual_simplex.hpp) — reports "
               "NUMERICAL_ERROR otherwise rather than a general dual phase 1\n"
               "  --solver managed   the \"concurrent solve manager\": tries dual, checker-\n"
               "                     verifies it, falls back to revised otherwise — see\n"
               "                     simplex/solve_manager.hpp for why this is a fallback\n"
               "                     chain, not real concurrency, for now\n"
               "  --presolve         fixed-variable + empty-column removal before solving "
               "(opt-in; only --solver revised uses it)\n";
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
  bool use_presolve = false;
  for (int i = 2; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--solution" && i + 1 < argc) {
      solution_path = argv[++i];
    } else if (arg == "--solver" && i + 1 < argc) {
      solver_name = argv[++i];
    } else if (arg == "--presolve") {
      use_presolve = true;
    }
  }
  if (solver_name != "revised" && solver_name != "dense" && solver_name != "dual" &&
      solver_name != "managed") {
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

  inferno::core::Solution solution;
  std::string winner;
  if (use_presolve && solver_name == "revised") {
    auto pre = inferno::presolve::Presolve(problem);
    if (pre.infeasible) {
      solution.status = inferno::core::SolveStatus::kInfeasible;
    } else {
      inferno::core::Solution reduced_solution = inferno::simplex::SolveRevised(pre.reduced);
      solution = inferno::presolve::Postsolve(problem, pre, reduced_solution);
    }
  } else if (solver_name == "dense") {
    solution = inferno::simplex::SolveDense(problem);
  } else if (solver_name == "dual") {
    solution = inferno::simplex::SolveDual(problem);
  } else if (solver_name == "managed") {
    auto managed = inferno::simplex::SolveManaged(problem);
    solution = managed.solution;
    winner = managed.winner;
  } else {
    solution = inferno::simplex::SolveRevised(problem);
  }

  auto end = std::chrono::steady_clock::now();
  double elapsed_s = std::chrono::duration<double>(end - start).count();

  std::cout << problem.name << " solver=" << solver_name << (use_presolve ? "+presolve" : "")
            << (winner.empty() ? "" : " winner=" + winner)
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
