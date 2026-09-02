// Phase 2.1 smoke test for the real revised simplex: same fixtures as
// tests/smoke_test.cpp's dense-simplex checks (so both solvers are held to
// the same bar), plus regression cases for bugs the dense simplex already
// had fixed for it (FX-variable complementarity) that a fresh
// implementation could easily reintroduce.

#include <cmath>
#include <iostream>
#include <string>

#include "checker/checker.hpp"
#include "core/lp_problem.hpp"
#include "core/sparse.hpp"
#include "simplex/revised_simplex.hpp"

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

using namespace inferno::core;
using inferno::simplex::SolveRevised;

// minimize -x - y  s.t.  x + y <= 4, 0 <= x <= 2, 0 <= y <= 3 -> optimum -4.
void TestHandBuiltLp() {
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

  Solution sol = SolveRevised(p);
  Expect(sol.status == SolveStatus::kOptimal, "revised: hand-built LP solves to optimal");
  if (sol.status == SolveStatus::kOptimal) {
    ExpectNear(sol.objective_value, -4.0, 1e-6, "revised: hand-built LP objective");
    auto check = inferno::checker::VerifySolution(p, sol);
    Expect(check.passed, "revised: hand-built LP checker passes: " + check.message);
  }
}

void TestUnboundedNoRows() {
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

  Solution sol = SolveRevised(p);
  Expect(sol.status == SolveStatus::kUnbounded, "revised: no-constraint LP detected as unbounded");
}

void TestFixedVariableComplementarityIsExempt() {
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

  Solution sol = SolveRevised(p);
  Expect(sol.status == SolveStatus::kOptimal, "revised: fixed-variable LP solves to optimal");
  if (sol.status == SolveStatus::kOptimal) {
    ExpectNear(sol.x[0], 5.0, 1e-9, "revised: fixed variable held at its bound");
    auto check = inferno::checker::VerifySolution(p, sol);
    Expect(check.passed, "revised: fixed-variable checker passes: " + check.message);
  }
}

// A range-constrained, multi-row LP with a mix of L/G/E rows and finite
// upper bounds — exercises RANGES-style row bounds and a real pivot
// sequence (not solvable in zero iterations from the slack basis).
//   minimize 2x + 3y - z
//   s.t.  1 <= x + y      <= 6      (range row)
//         x - y + 2z       = 3      (equality row)
//                y +  z   <= 5
//         0 <= x <= 4, 0 <= y <= 4, 0 <= z <= 4
void TestMultiRowRangeLp() {
  LpProblem p;
  p.name = "multi_row";
  p.num_rows = 3;
  p.num_cols = 3;

  CscBuilder builder(3, 3);
  builder.AddEntry(0, 0, 1.0);  // x in row0
  builder.AddEntry(0, 1, 1.0);  // x in row1
  builder.AddEntry(1, 0, 1.0);  // y in row0
  builder.AddEntry(1, 1, -1.0); // y in row1
  builder.AddEntry(1, 2, 1.0);  // y in row2
  builder.AddEntry(2, 1, 2.0);  // z in row1
  builder.AddEntry(2, 2, 1.0);  // z in row2
  p.a = std::move(builder).Build();

  p.row_lo = {1.0, 3.0, -kInfinity};
  p.row_hi = {6.0, 3.0, 5.0};
  p.row_names = {"R0", "R1", "R2"};
  p.col_lo = {0.0, 0.0, 0.0};
  p.col_hi = {4.0, 4.0, 4.0};
  p.obj = {2.0, 3.0, -1.0};
  p.col_names = {"x", "y", "z"};

  Solution sol = SolveRevised(p);
  Expect(sol.status == SolveStatus::kOptimal, "revised: multi-row range LP solves to optimal");
  if (sol.status == SolveStatus::kOptimal) {
    auto check = inferno::checker::VerifySolution(p, sol);
    Expect(check.passed, "revised: multi-row range LP checker passes: " + check.message);
  }
}

}  // namespace

int main() {
  TestHandBuiltLp();
  TestUnboundedNoRows();
  TestFixedVariableComplementarityIsExempt();
  TestMultiRowRangeLp();

  if (g_failures > 0) {
    std::cerr << g_failures << " test(s) failed\n";
    return 1;
  }
  std::cout << "all revised simplex tests passed\n";
  return 0;
}
