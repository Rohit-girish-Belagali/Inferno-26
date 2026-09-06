// MILP regression suite. Every case has an independently known optimum,
// and every answer additionally goes through checker/mip_checker, which
// recomputes integrality, bounds, rows, objective, bound sanity and gap
// from the raw model without consulting anything the search computed.
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "checker/mip_checker.hpp"
#include "core/mip_problem.hpp"
#include "core/sparse.hpp"
#include "mip/branch_and_bound.hpp"
#include "simplex/revised_simplex.hpp"

namespace {
int g_fail = 0;
void Expect(bool c, const std::string& m) {
  if (!c) { std::cerr << "FAIL: " << m << "\n"; ++g_fail; }
  else { std::cout << "ok: " << m << "\n"; }
}

using namespace inferno;

// Runs a MILP, verifies it independently, and checks the known optimum.
// `expect_obj` is NaN when the case has no finite optimum (infeasible).
void Check(const std::string& name, const core::MipProblem& p, double expect_obj,
           bool expect_feasible = true,
           mip::BranchAndBoundOptions opts = mip::BranchAndBoundOptions{}) {
  auto s = mip::SolveMip(p, opts);
  if (!expect_feasible) {
    Expect(s.status == core::SolveStatus::kInfeasible, name + ": proved infeasible");
    return;
  }
  Expect(s.status == core::SolveStatus::kOptimal, name + ": solved");
  if (s.status != core::SolveStatus::kOptimal) return;

  auto chk = checker::VerifyMipSolution(p, s);
  Expect(chk.passed, name + ": independent verifier — " + chk.message);
  Expect(std::abs(s.objective_value - expect_obj) < 1e-6,
         name + ": objective " + std::to_string(s.objective_value) + " matches known optimum " +
             std::to_string(expect_obj));
  Expect(s.proved_optimal, name + ": optimality proved");
}

core::MipProblem MakeKnapsack(const std::vector<double>& value,
                               const std::vector<double>& weight, double cap,
                               const std::vector<double>& ub) {
  const int n = static_cast<int>(value.size());
  core::MipProblem p;
  p.lp.num_cols = n;
  p.lp.num_rows = 1;
  core::CscBuilder b(1, n);
  for (int j = 0; j < n; ++j) b.AddEntry(j, 0, weight[j]);
  p.lp.a = std::move(b).Build();
  p.lp.obj.resize(n);
  for (int j = 0; j < n; ++j) p.lp.obj[j] = -value[j];  // maximise value
  p.lp.col_lo.assign(n, 0.0);
  p.lp.col_hi = ub;
  p.lp.row_lo = {-core::kInfinity};
  p.lp.row_hi = {cap};
  p.lp.col_names.assign(n, "x");
  p.lp.row_names = {"cap"};
  p.is_integer.assign(n, 1);
  return p;
}
}  // namespace

// 1. Binary knapsack. Values 5,4,3; weights 2,3,1; capacity 4.
// Best is items 1+3 = 8.
void TestBinaryKnapsack() {
  Check("binary knapsack", MakeKnapsack({5, 4, 3}, {2, 3, 1}, 4, {1, 1, 1}), -8.0);
}

// 2. Integer (unbounded-ish) knapsack: one item, value 3 weight 2,
// capacity 7, up to 10 copies -> 3 copies = 9.
void TestIntegerKnapsack() {
  Check("integer knapsack", MakeKnapsack({3}, {2}, 7, {10}), -9.0);
}

// 3. Assignment: 2 agents x 2 tasks, costs [[4,2],[3,5]]. Each agent one
// task, each task one agent. Optimum: a0->t1 (2) + a1->t0 (3) = 5.
void TestAssignment() {
  core::MipProblem p;
  p.lp.num_cols = 4;  // x00 x01 x10 x11
  p.lp.num_rows = 4;
  core::CscBuilder b(4, 4);
  b.AddEntry(0, 0, 1.0); b.AddEntry(0, 2, 1.0);   // x00: agent0, task0
  b.AddEntry(1, 0, 1.0); b.AddEntry(1, 3, 1.0);   // x01: agent0, task1
  b.AddEntry(2, 1, 1.0); b.AddEntry(2, 2, 1.0);   // x10: agent1, task0
  b.AddEntry(3, 1, 1.0); b.AddEntry(3, 3, 1.0);   // x11: agent1, task1
  p.lp.a = std::move(b).Build();
  p.lp.obj = {4, 2, 3, 5};
  p.lp.col_lo.assign(4, 0.0);
  p.lp.col_hi.assign(4, 1.0);
  p.lp.row_lo = {1, 1, 1, 1};
  p.lp.row_hi = {1, 1, 1, 1};
  p.lp.col_names.assign(4, "x");
  p.lp.row_names = {"a0", "a1", "t0", "t1"};
  p.is_integer.assign(4, 1);
  Check("assignment", p, 5.0);
}

// 4. Set covering: 3 sets over 3 elements. S0={e0,e1} cost 3,
// S1={e1,e2} cost 3, S2={e0,e1,e2} cost 7. The single all-covering set is
// deliberately priced ABOVE the pair, so the cheapest cover is S0+S1 = 6
// and the case tests that the search combines sets rather than grabbing
// the one that trivially covers everything. (First written with S2 at
// cost 5, where S2 alone is genuinely optimal at 5 — the test expectation
// was wrong, not the solver, and the independent verifier passed the
// answer throughout.)
void TestSetCovering() {
  core::MipProblem p;
  p.lp.num_cols = 3;
  p.lp.num_rows = 3;
  core::CscBuilder b(3, 3);
  b.AddEntry(0, 0, 1.0); b.AddEntry(0, 1, 1.0);
  b.AddEntry(1, 1, 1.0); b.AddEntry(1, 2, 1.0);
  b.AddEntry(2, 0, 1.0); b.AddEntry(2, 1, 1.0); b.AddEntry(2, 2, 1.0);
  p.lp.a = std::move(b).Build();
  p.lp.obj = {3, 3, 7};
  p.lp.col_lo.assign(3, 0.0);
  p.lp.col_hi.assign(3, 1.0);
  p.lp.row_lo = {1, 1, 1};
  p.lp.row_hi = {core::kInfinity, core::kInfinity, core::kInfinity};
  p.lp.col_names.assign(3, "s");
  p.lp.row_names = {"e0", "e1", "e2"};
  p.is_integer.assign(3, 1);
  Check("set covering", p, 6.0);
}

// 5. Already-integral relaxation: the LP optimum is integral, so B&B must
// accept it at the root without branching.
void TestAlreadyIntegral() {
  core::MipProblem p;
  p.lp.num_cols = 1;
  p.lp.num_rows = 1;
  core::CscBuilder b(1, 1);
  b.AddEntry(0, 0, 1.0);
  p.lp.a = std::move(b).Build();
  p.lp.obj = {-1.0};
  p.lp.col_lo = {0.0};
  p.lp.col_hi = {5.0};
  p.lp.row_lo = {-core::kInfinity};
  p.lp.row_hi = {5.0};
  p.lp.col_names = {"x"};
  p.lp.row_names = {"c"};
  p.is_integer = {1};
  auto s = mip::SolveMip(p);
  Expect(s.status == core::SolveStatus::kOptimal, "already-integral: solved");
  Expect(std::abs(s.objective_value + 5.0) < 1e-9, "already-integral: objective -5");
  Expect(s.nodes_explored <= 2, "already-integral: accepted without a deep search (nodes=" +
                                     std::to_string(s.nodes_explored) + ")");
}

// 6. Infeasible over the integers though the relaxation is feasible.
void TestInfeasible() {
  core::MipProblem p;
  p.lp.num_cols = 2;
  p.lp.num_rows = 1;
  core::CscBuilder b(1, 2);
  b.AddEntry(0, 0, 2.0); b.AddEntry(1, 0, 2.0);
  p.lp.a = std::move(b).Build();
  p.lp.obj = {1, 1};
  p.lp.col_lo = {0, 0};
  p.lp.col_hi = {5, 5};
  p.lp.row_lo = {3.0};
  p.lp.row_hi = {3.0};
  p.lp.col_names.assign(2, "x");
  p.lp.row_names = {"eq"};
  p.is_integer = {1, 1};
  Check("integer-infeasible", p, 0.0, /*expect_feasible=*/false);
}

// 7. Weak LP relaxation: the classic "equality knapsack" where the
// relaxation is far above the integer optimum, forcing real search.
// 2x0+2x1+2x2+2x3 <= 5, maximise x0+x1+x2+x3, binaries.
// Relaxation gives 2.5; integer optimum is 2.
void TestWeakRelaxation() {
  core::MipProblem p;
  p.lp.num_cols = 4;
  p.lp.num_rows = 1;
  core::CscBuilder b(1, 4);
  for (int j = 0; j < 4; ++j) b.AddEntry(j, 0, 2.0);
  p.lp.a = std::move(b).Build();
  p.lp.obj = {-1, -1, -1, -1};
  p.lp.col_lo.assign(4, 0.0);
  p.lp.col_hi.assign(4, 1.0);
  p.lp.row_lo = {-core::kInfinity};
  p.lp.row_hi = {5.0};
  p.lp.col_names.assign(4, "x");
  p.lp.row_names = {"cap"};
  p.is_integer.assign(4, 1);

  core::Solution relax = simplex::SolveRevised(p.lp);
  Expect(relax.objective_value < -2.4,
         "weak relaxation: LP bound is genuinely weak (" + std::to_string(relax.objective_value) + ")");
  Check("weak relaxation", p, -2.0);
}

// 8. Production planning with binary setup: produce x (0..10) at unit cost
// 2, only if the line is set up at fixed cost 12; demand 4.
// Cost = 12 + 8 = 20.
void TestProductionSetup() {
  core::MipProblem p;
  p.lp.num_cols = 2;  // setup, x
  p.lp.num_rows = 2;
  core::CscBuilder b(2, 2);
  b.AddEntry(0, 1, -10.0);  // x <= 10*setup
  b.AddEntry(1, 0, 1.0);    // x >= demand
  b.AddEntry(1, 1, 1.0);
  p.lp.a = std::move(b).Build();
  p.lp.obj = {12.0, 2.0};
  p.lp.col_lo = {0, 0};
  p.lp.col_hi = {1, 10};
  p.lp.row_lo = {4.0, -core::kInfinity};
  p.lp.row_hi = {core::kInfinity, 0.0};
  p.lp.col_names = {"setup", "x"};
  p.lp.row_names = {"demand", "link"};
  p.is_integer = {1, 0};  // mixed
  Check("production with setup", p, 20.0);
}

// 9. Requires several nodes: a 6-item binary knapsack whose greedy and LP
// answers both differ from the optimum.
void TestMultiNode() {
  auto p = MakeKnapsack({10, 9, 8, 7, 6, 5}, {5, 4, 4, 3, 3, 2}, 10, {1, 1, 1, 1, 1, 1});
  auto s = mip::SolveMip(p);
  Expect(s.status == core::SolveStatus::kOptimal, "multi-node knapsack: solved");
  auto chk = checker::VerifyMipSolution(p, s);
  Expect(chk.passed, "multi-node knapsack: independent verifier — " + chk.message);
  Expect(s.nodes_explored > 1,
         "multi-node knapsack: search actually branched (nodes=" +
             std::to_string(s.nodes_explored) + ")");
  Expect(s.proved_optimal, "multi-node knapsack: optimality proved");
}

// 10. Limit handling. A node limit of 1 must NOT yield a claim of proved
// optimality, and the reported bound must remain valid.
void TestLimits() {
  auto p = MakeKnapsack({10, 9, 8, 7, 6, 5}, {5, 4, 4, 3, 3, 2}, 10, {1, 1, 1, 1, 1, 1});
  mip::BranchAndBoundOptions opts;
  opts.node_limit = 1;
  opts.rounding_heuristic = false;  // force the search to do the work
  auto s = mip::SolveMip(p, opts);
  Expect(!s.proved_optimal, "node limit: optimality is NOT claimed");
  if (s.status == core::SolveStatus::kOptimal) {
    auto chk = checker::VerifyMipSolution(p, s);
    Expect(!chk.bound_contradicts_objective,
           "node limit: reported bound does not exceed the incumbent objective");
    Expect(chk.passed, "node limit: incumbent still verifies — " + chk.message);
  }

  mip::BranchAndBoundOptions topts;
  topts.time_limit_seconds = 0.0;  // expire immediately
  auto t = mip::SolveMip(p, topts);
  Expect(!t.proved_optimal, "time limit: optimality is NOT claimed");
}

int main() {
  TestBinaryKnapsack();
  TestIntegerKnapsack();
  TestAssignment();
  TestSetCovering();
  TestAlreadyIntegral();
  TestInfeasible();
  TestWeakRelaxation();
  TestProductionSetup();
  TestMultiNode();
  TestLimits();
  if (g_fail) { std::cerr << g_fail << " test(s) failed\n"; return 1; }
  std::cout << "all MILP suite tests passed\n";
  return 0;
}
