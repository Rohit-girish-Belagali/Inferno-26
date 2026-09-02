// Phase 2.2 presolve test. BUILD_PLAN_V2.md's own standing rule: "Presolve
// on/off equivalence is an automated test, not a habit" — this is that
// test. It solves real Netlib instances both with and without presolve and
// checks the objective matches and the independent checker passes both
// ways, plus a couple of hand-built cases for the specific reduction
// kinds and for presolve-detected infeasibility.

#include <cmath>
#include <iostream>
#include <string>

#include "checker/checker.hpp"
#include "core/lp_problem.hpp"
#include "core/sparse.hpp"
#include "io/mps_reader.hpp"
#include "presolve/presolve.hpp"
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
using inferno::presolve::Postsolve;
using inferno::presolve::Presolve;
using inferno::simplex::SolveRevised;

// minimize -x - 2z  s.t.  x + z <= 4,  0 <= x <= 2,  z fixed at 1,
// y present but empty (appears in no row), cost +3 (favors y at its lower
// bound, 0), so the true optimum is x=2 (bound), z=1 (fixed), y=0.
// Objective = -2 - 2 = -4.
void TestFixedAndEmptyColumn() {
  LpProblem p;
  p.name = "fixed_and_empty";
  p.num_rows = 1;
  p.num_cols = 3;  // x, z, y

  CscBuilder builder(1, 3);
  builder.AddEntry(0, 0, 1.0);  // x in row 0
  builder.AddEntry(1, 0, 1.0);  // z in row 0
  // y (col 2) has no entries: an empty column.
  p.a = std::move(builder).Build();

  p.row_lo = {-kInfinity};
  p.row_hi = {4.0};
  p.row_names = {"C1"};
  p.col_lo = {0.0, 1.0, 0.0};
  p.col_hi = {2.0, 1.0, kInfinity};
  p.obj = {-1.0, -2.0, 3.0};
  p.col_names = {"x", "z", "y"};

  auto pre = Presolve(p);
  Expect(!pre.infeasible, "fixed+empty-column case: presolve does not report infeasible");
  Expect(pre.reduced.num_cols == 1, "fixed+empty-column case: only x survives presolve (got " +
                                         std::to_string(pre.reduced.num_cols) + " columns)");

  Solution reduced_sol = SolveRevised(pre.reduced);
  Expect(reduced_sol.status == SolveStatus::kOptimal, "reduced problem solves to optimal");

  Solution sol = Postsolve(p, pre, reduced_sol);
  Expect(sol.status == SolveStatus::kOptimal, "postsolved solution reports optimal");
  if (sol.status == SolveStatus::kOptimal) {
    ExpectNear(sol.objective_value, -4.0, 1e-9, "postsolved objective");
    ExpectNear(sol.x[1], 1.0, 1e-9, "fixed variable z restored to its fixed value");
    ExpectNear(sol.x[2], 0.0, 1e-9, "empty column y restored to its favorable bound");
    auto check = inferno::checker::VerifySolution(p, sol);
    Expect(check.passed, "postsolved solution passes the independent checker: " + check.message);
  }
}

// x fixed at 5, but the only row forces x <= 3 — infeasible once x's fixed
// value is folded in, and presolve should catch it without ever calling
// the solver.
void TestPresolveDetectsInfeasibility() {
  LpProblem p;
  p.name = "infeasible_after_fixing";
  p.num_rows = 1;
  p.num_cols = 1;

  CscBuilder builder(1, 1);
  builder.AddEntry(0, 0, 1.0);
  p.a = std::move(builder).Build();

  p.row_lo = {-kInfinity};
  p.row_hi = {3.0};
  p.row_names = {"C1"};
  p.col_lo = {5.0};
  p.col_hi = {5.0};
  p.obj = {1.0};
  p.col_names = {"x"};

  auto pre = Presolve(p);
  Expect(pre.infeasible, "fixing x=5 into a row capping x<=3 is caught as infeasible by presolve");
}

// The plan's own standing rule: presolve on/off equivalence, checked
// against real data — several Netlib instances are known (from earlier
// debugging) to actually contain FX-bound (fixed) variables, so this
// exercises the real reduction path, not a no-op.
void TestNetlibEquivalence(const std::string& name) {
  std::string path = "bench/netlib/mps/" + name + ".mps";
  LpProblem p;
  try {
    p = inferno::io::ReadMps(path);
  } catch (const std::exception& e) {
    std::cout << "skip (no local Netlib set): " << name << " — " << e.what() << "\n";
    return;
  }

  Solution direct = SolveRevised(p);
  auto pre = Presolve(p);
  Expect(!pre.infeasible, name + ": presolve does not falsely report infeasible");
  if (pre.infeasible) return;

  Solution reduced_sol = SolveRevised(pre.reduced);
  Solution postsolved = Postsolve(p, pre, reduced_sol);

  Expect(direct.status == postsolved.status,
         name + ": presolve on/off status matches (direct=" + std::to_string((int)direct.status) +
             " postsolved=" + std::to_string((int)postsolved.status) + ")");
  if (direct.status == SolveStatus::kOptimal && postsolved.status == SolveStatus::kOptimal) {
    ExpectNear(postsolved.objective_value, direct.objective_value,
               1e-6 * std::max(1.0, std::abs(direct.objective_value)),
               name + ": presolve on/off objective matches");
    auto check = inferno::checker::VerifySolution(p, postsolved);
    Expect(check.passed, name + ": postsolved solution passes the independent checker: " +
                              check.message);
  }
}

}  // namespace

int main() {
  TestFixedAndEmptyColumn();
  TestPresolveDetectsInfeasibility();
  // capri, finnis, recipe, shell, standata, standgub, standmps and
  // vtp.base are known (from earlier checker debugging this session) to
  // contain real FX-bound variables, so presolve actually does something
  // on these rather than a vacuous no-op.
  for (const char* name : {"capri", "finnis", "recipe", "shell", "standata", "standgub",
                            "standmps", "vtp.base", "afiro", "blend"}) {
    TestNetlibEquivalence(name);
  }

  if (g_failures > 0) {
    std::cerr << g_failures << " test(s) failed\n";
    return 1;
  }
  std::cout << "all presolve tests passed\n";
  return 0;
}
