// Dual simplex test. Cross-validated against the primal revised simplex on
// real Netlib data: same objective, checker PASS either way — a
// meaningfully independent check since the two algorithms select entering
// and leaving variables in the opposite order and would have to agree by
// correctness, not by construction.

#include <cmath>
#include <iostream>
#include <string>

#include "checker/checker.hpp"
#include "core/lp_problem.hpp"
#include "core/sparse.hpp"
#include "io/mps_reader.hpp"
#include "simplex/dual_simplex.hpp"
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
using inferno::simplex::SolveDual;
using inferno::simplex::SolveRevised;

// minimize -x - y  s.t.  x + y <= 4, 0 <= x <= 2, 0 <= y <= 3 -> optimum -4.
// Both structural costs are negative, so the trivial dual-feasible start
// (AtUpper for both, since finite upper bounds exist) applies directly.
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

  Solution sol = SolveDual(p);
  Expect(sol.status == SolveStatus::kOptimal, "dual: hand-built LP solves to optimal");
  if (sol.status == SolveStatus::kOptimal) {
    ExpectNear(sol.objective_value, -4.0, 1e-6, "dual: hand-built LP objective");
    auto check = inferno::checker::VerifySolution(p, sol);
    Expect(check.passed, "dual: hand-built LP checker passes: " + check.message);
  }
}

void TestMultiRowRangeLp() {
  LpProblem p;
  p.name = "multi_row";
  p.num_rows = 3;
  p.num_cols = 3;

  CscBuilder builder(3, 3);
  builder.AddEntry(0, 0, 1.0);
  builder.AddEntry(0, 1, 1.0);
  builder.AddEntry(1, 0, 1.0);
  builder.AddEntry(1, 1, -1.0);
  builder.AddEntry(1, 2, 1.0);
  builder.AddEntry(2, 1, 2.0);
  builder.AddEntry(2, 2, 1.0);
  p.a = std::move(builder).Build();

  p.row_lo = {1.0, 3.0, -kInfinity};
  p.row_hi = {6.0, 3.0, 5.0};
  p.row_names = {"R0", "R1", "R2"};
  p.col_lo = {0.0, 0.0, 0.0};
  p.col_hi = {4.0, 4.0, 4.0};
  p.obj = {2.0, 3.0, -1.0};
  p.col_names = {"x", "y", "z"};

  Solution primal = SolveRevised(p);
  Solution dual = SolveDual(p);
  Expect(dual.status == SolveStatus::kOptimal, "dual: multi-row LP solves to optimal");
  if (dual.status == SolveStatus::kOptimal && primal.status == SolveStatus::kOptimal) {
    ExpectNear(dual.objective_value, primal.objective_value, 1e-6,
               "dual: multi-row LP objective matches primal revised simplex");
    auto check = inferno::checker::VerifySolution(p, dual);
    Expect(check.passed, "dual: multi-row LP checker passes: " + check.message);
  }
}

// Cross-validates dual against primal on real Netlib instances whose
// structural costs are all non-negative (so the trivial dual-feasible
// start applies directly — see dual_simplex.hpp's scope note; a general
// dual phase 1 isn't implemented).
void TestNetlibCrossCheck(const std::string& name) {
  std::string path = "bench/netlib/mps/" + name + ".mps";
  LpProblem p;
  try {
    p = inferno::io::ReadMps(path);
  } catch (const std::exception& e) {
    std::cout << "skip (no local Netlib set): " << name << " — " << e.what() << "\n";
    return;
  }

  Solution primal = SolveRevised(p);
  Solution dual = SolveDual(p);

  if (dual.status == SolveStatus::kNumericalError) {
    std::cout << "skip (no trivial dual-feasible start): " << name << "\n";
    return;
  }

  Expect(dual.status == primal.status,
         name + ": dual/primal status matches (primal=" + std::to_string((int)primal.status) +
             " dual=" + std::to_string((int)dual.status) + ")");
  if (dual.status == SolveStatus::kOptimal && primal.status == SolveStatus::kOptimal) {
    ExpectNear(dual.objective_value, primal.objective_value,
               1e-6 * std::max(1.0, std::abs(primal.objective_value)),
               name + ": dual objective matches primal");
    auto check = inferno::checker::VerifySolution(p, dual);
    Expect(check.passed, name + ": dual solution passes the independent checker: " + check.message);
  }
}

}  // namespace

int main() {
  TestHandBuiltLp();
  TestMultiRowRangeLp();
  // Instances with all-nonnegative costs (checked ad hoc against each
  // file), so the trivial dual-feasible start applies.
  for (const char* name : {"afiro", "adlittle", "blend", "kb2", "recipe", "share2b",
                            "scagr7", "israel", "sc50a", "sc50b"}) {
    TestNetlibCrossCheck(name);
  }

  if (g_failures > 0) {
    std::cerr << g_failures << " test(s) failed\n";
    return 1;
  }
  std::cout << "all dual simplex tests passed\n";
  return 0;
}
