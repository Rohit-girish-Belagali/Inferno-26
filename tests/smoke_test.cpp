// Phase 1.1 smoke test: exercises the MPS reader, dense proof-of-life
// simplex and independent checker together, plus a couple of edge cases
// (unbounded, no-constraint problems) that are cheap to get wrong.
//
// Deliberately dependency-free (no GoogleTest) — Phase 1.1 has no test
// framework decision made yet, and pulling one in is out of scope here.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include "checker/checker.hpp"
#include "core/lp_problem.hpp"
#include "core/sparse.hpp"
#include "io/mps_reader.hpp"
#include "simplex/dense_simplex.hpp"

namespace {

int g_failures = 0;

void Expect(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << "\n";
    ++g_failures;
  } else {
    std::cout << "ok: " << message << "\n";
  }
}

void ExpectNear(double actual, double expected, double tol, const std::string& message) {
  Expect(std::abs(actual - expected) <= tol,
         message + " (expected " + std::to_string(expected) + ", got " + std::to_string(actual) +
             ")");
}

// minimize -x - y  s.t.  x + y <= 4, 0 <= x <= 2, 0 <= y <= 3
// Optimal objective is -4 (alternate optima exist on the x+y=4 edge, so we
// only check objective value and feasibility, not a specific vertex).
void TestHandBuiltLp() {
  using namespace inferno::core;

  LpProblem p;
  p.name = "hand_built";
  p.num_rows = 1;
  p.num_cols = 2;

  CscBuilder builder(1, 2);
  builder.AddEntry(0, 0, 1.0);
  builder.AddEntry(1, 0, 1.0);
  p.a = std::move(builder).Build();

  p.row_lo = {-kInfinity};
  p.row_hi = {4.0};
  p.row_names = {"C1"};

  p.col_lo = {0.0, 0.0};
  p.col_hi = {2.0, 3.0};
  p.obj = {-1.0, -1.0};
  p.col_names = {"X", "Y"};

  Solution sol = inferno::simplex::SolveDense(p);
  Expect(sol.status == SolveStatus::kOptimal, "hand-built LP solves to optimal");
  if (sol.status == SolveStatus::kOptimal) {
    ExpectNear(sol.objective_value, -4.0, 1e-6, "hand-built LP objective");
    auto check = inferno::checker::VerifySolution(p, sol);
    Expect(check.passed, "hand-built LP checker passes: " + check.message);
  }
}

// A fixed variable (lo == hi) is "at both bounds" simultaneously and is
// exempt from the complementarity condition — its reduced cost may be any
// sign at optimality. Regression test for a checker bug where FX-bound
// variables were flagged as complementarity violations regardless of the
// sign of their (perfectly valid) nonzero reduced cost.
void TestFixedVariableComplementarityIsExempt() {
  using namespace inferno::core;

  LpProblem p;
  p.name = "fixed_var";
  p.num_rows = 0;
  p.num_cols = 1;

  CscMatrix a;
  a.rows = 0;
  a.cols = 1;
  a.col_ptr = {0, 0};
  p.a = a;

  p.col_lo = {5.0};
  p.col_hi = {5.0};
  p.obj = {7.0};
  p.col_names = {"Z"};

  Solution sol;
  sol.status = SolveStatus::kOptimal;
  sol.x = {5.0};
  sol.row_activity = {};
  sol.y = {};
  sol.reduced_cost = {7.0};  // == obj[0] - A^T y with no rows: a valid, nonzero z.
  sol.objective_value = 35.0;

  auto check = inferno::checker::VerifySolution(p, sol);
  Expect(check.passed, "fixed variable with nonzero reduced cost passes checker: " + check.message);
}

// minimize -x, x >= 0, no upper bound, no constraints at all (num_rows =
// 0). Exercises the m = 0 edge case in the dense simplex.
void TestUnboundedNoRows() {
  using namespace inferno::core;

  LpProblem p;
  p.name = "unbounded";
  p.num_rows = 0;
  p.num_cols = 1;

  CscMatrix a;
  a.rows = 0;
  a.cols = 1;
  a.col_ptr = {0, 0};
  p.a = a;

  p.col_lo = {0.0};
  p.col_hi = {kInfinity};
  p.obj = {-1.0};
  p.col_names = {"X"};

  Solution sol = inferno::simplex::SolveDense(p);
  Expect(sol.status == SolveStatus::kUnbounded, "no-constraint LP detected as unbounded");
}

// Round-trips a tiny hand-written MPS file through the reader, then solves
// it and checks against the same known optimum as TestHandBuiltLp — this
// is the closest thing Phase 1.1 has to Checkpoint 1.1's afiro check
// without shipping the Netlib set into the repo.
void TestMpsRoundTrip() {
  const std::string path = "tiny_test.mps";
  {
    std::ofstream out(path);
    out << "NAME TINY\n"
           "ROWS\n"
           " N COST\n"
           " L C1\n"
           "COLUMNS\n"
           " X COST -1.0 C1 1.0\n"
           " Y COST -1.0 C1 1.0\n"
           "RHS\n"
           " RHS C1 4.0\n"
           "BOUNDS\n"
           " UP BND X 2.0\n"
           " UP BND Y 3.0\n"
           "ENDATA\n";
  }

  inferno::core::LpProblem p = inferno::io::ReadMps(path);
  Expect(p.num_rows == 1 && p.num_cols == 2, "MPS round-trip parses expected shape");

  inferno::core::Solution sol = inferno::simplex::SolveDense(p);
  Expect(sol.status == inferno::core::SolveStatus::kOptimal, "MPS round-trip solves to optimal");
  if (sol.status == inferno::core::SolveStatus::kOptimal) {
    ExpectNear(sol.objective_value, -4.0, 1e-6, "MPS round-trip objective");
    auto check = inferno::checker::VerifySolution(p, sol);
    Expect(check.passed, "MPS round-trip checker passes: " + check.message);
  }

  std::remove(path.c_str());
}

}  // namespace

int main() {
  TestHandBuiltLp();
  TestFixedVariableComplementarityIsExempt();
  TestUnboundedNoRows();
  TestMpsRoundTrip();

  if (g_failures > 0) {
    std::cerr << g_failures << " test(s) failed\n";
    return 1;
  }
  std::cout << "all tests passed\n";
  return 0;
}
