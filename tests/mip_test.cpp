// Branch-and-bound tests. Each case has an optimum derivable by hand, and
// each is chosen so the LP relaxation is FRACTIONAL — otherwise the test
// would pass without branching ever happening and would prove nothing
// about the search.
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "core/mip_problem.hpp"
#include "core/sparse.hpp"
#include "mip/branch_and_bound.hpp"
#include "simplex/revised_simplex.hpp"

namespace {
int g_failures = 0;
void Expect(bool c, const std::string& m) {
  if (!c) { std::cerr << "FAIL: " << m << "\n"; ++g_failures; }
  else { std::cout << "ok: " << m << "\n"; }
}
void ExpectNear(double a, double e, double t, const std::string& m) {
  Expect(std::abs(a - e) <= t,
         m + " (expected " + std::to_string(e) + ", got " + std::to_string(a) + ")");
}
}  // namespace

using namespace inferno;

// 0/1 knapsack: maximise 5a + 4b + 3c s.t. 2a + 3b + c <= 4, all binary.
// LP relaxation is fractional. Best integer point is a=1, c=1 -> 8
// (a and b together need capacity 5 > 4).
void TestKnapsack() {
  core::MipProblem p;
  p.lp.num_cols = 3;
  p.lp.num_rows = 1;
  core::CscBuilder b(1, 3);
  b.AddEntry(0, 0, 2.0); b.AddEntry(1, 0, 3.0); b.AddEntry(2, 0, 1.0);
  p.lp.a = std::move(b).Build();
  p.lp.obj = {-5.0, -4.0, -3.0};
  p.lp.col_lo = {0, 0, 0};
  p.lp.col_hi = {1, 1, 1};
  p.lp.row_lo = {-core::kInfinity};
  p.lp.row_hi = {4.0};
  p.lp.col_names = {"a", "b", "c"};
  p.lp.row_names = {"cap"};
  p.is_integer = {1, 1, 1};

  // Confirm the relaxation really is fractional, so branching is exercised.
  core::Solution relax = simplex::SolveRevised(p.lp);
  bool frac = false;
  for (double v : relax.x) if (std::abs(v - std::round(v)) > 1e-6) frac = true;
  Expect(frac, "knapsack LP relaxation is fractional (so branching is required)");

  auto s = mip::SolveMip(p);
  Expect(s.status == core::SolveStatus::kOptimal, "knapsack solves");
  ExpectNear(s.objective_value, -8.0, 1e-6, "knapsack optimal objective is -8");
  Expect(s.proved_optimal, "knapsack optimality is proved, not merely reached");
  for (int j = 0; j < 3; ++j) {
    Expect(std::abs(s.x[j] - std::round(s.x[j])) < 1e-6, "knapsack x" + std::to_string(j) + " is integral");
  }
  Expect(s.gap < 1e-6, "knapsack gap closed");
}

// Facility/plant selection in miniature: two plants with fixed opening
// costs and capacities, demand 15. Opening only the big plant costs
// 100 + 15*2 = 130; only the small one cannot meet demand (cap 10); both
// costs 100 + 60 + 15*2ish. Optimum is the big plant alone at 130.
void TestPlantSelection() {
  // cols: open_big, open_small, prod_big, prod_small
  core::MipProblem p;
  p.lp.num_cols = 4;
  p.lp.num_rows = 3;
  core::CscBuilder b(3, 4);
  // row0: prod_big + prod_small >= 15  (demand)
  // row1: prod_big - 20*open_big <= 0  (capacity links production to opening)
  // row2: prod_small - 10*open_small <= 0
  b.AddEntry(0, 1, -20.0);
  b.AddEntry(1, 2, -10.0);
  b.AddEntry(2, 0, 1.0); b.AddEntry(2, 1, 1.0);
  b.AddEntry(3, 0, 1.0); b.AddEntry(3, 2, 1.0);
  p.lp.a = std::move(b).Build();
  p.lp.obj = {100.0, 60.0, 2.0, 2.0};
  p.lp.col_lo = {0, 0, 0, 0};
  p.lp.col_hi = {1, 1, 20, 10};
  p.lp.row_lo = {15.0, -core::kInfinity, -core::kInfinity};
  p.lp.row_hi = {core::kInfinity, 0.0, 0.0};
  p.lp.col_names = {"open_big", "open_small", "prod_big", "prod_small"};
  p.lp.row_names = {"demand", "cap_big", "cap_small"};
  p.is_integer = {1, 1, 0, 0};  // mixed: opening binary, production continuous

  auto s = mip::SolveMip(p);
  Expect(s.status == core::SolveStatus::kOptimal, "plant selection solves");
  ExpectNear(s.objective_value, 130.0, 1e-6, "plant selection optimal cost is 130");
  Expect(std::abs(s.x[0] - 1.0) < 1e-6, "the big plant is opened");
  Expect(std::abs(s.x[1]) < 1e-6, "the small plant is not opened");
  Expect(s.proved_optimal, "plant selection optimality is proved");
  // The mixed part: production is continuous and need not be integral.
  ExpectNear(s.x[2], 15.0, 1e-6, "big plant produces exactly the demand");
}

// An equality forcing a fractional relaxation: 2x + 2y = 3 has no integer
// solution at all, so the search must PROVE infeasibility rather than
// return the relaxation's fractional point.
void TestIntegerInfeasible() {
  core::MipProblem p;
  p.lp.num_cols = 2;
  p.lp.num_rows = 1;
  core::CscBuilder b(1, 2);
  b.AddEntry(0, 0, 2.0); b.AddEntry(1, 0, 2.0);
  p.lp.a = std::move(b).Build();
  p.lp.obj = {1.0, 1.0};
  p.lp.col_lo = {0, 0};
  p.lp.col_hi = {5, 5};
  p.lp.row_lo = {3.0};
  p.lp.row_hi = {3.0};
  p.lp.col_names = {"x", "y"};
  p.lp.row_names = {"eq"};
  p.is_integer = {1, 1};

  auto s = mip::SolveMip(p);
  Expect(s.status == core::SolveStatus::kInfeasible,
         "2x + 2y = 3 with x,y integer is proved infeasible, not answered fractionally");
}

int main() {
  TestKnapsack();
  TestPlantSelection();
  TestIntegerInfeasible();
  if (g_failures) { std::cerr << g_failures << " test(s) failed\n"; return 1; }
  std::cout << "all MIP tests passed\n";
  return 0;
}
