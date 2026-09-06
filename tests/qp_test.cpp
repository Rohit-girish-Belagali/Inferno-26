// Phase 4Q QP test. Each case has an optimum derivable by hand, so the
// test checks the ANSWER, not merely that ADMM terminated.
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "core/qp_problem.hpp"
#include "core/sparse.hpp"
#include "qp/admm.hpp"

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

// min 1/2(x0^2 + x1^2) - x0 - x1,  s.t. 0 <= x <= 10 (as rows).
// Unconstrained optimum is x = (1, 1), interior to the box, so the bounds
// are inactive and both duals must come out zero. Objective = -1.
void TestUnconstrainedInterior() {
  core::QpProblem p;
  p.name = "interior";
  p.num_cols = 2;
  p.num_rows = 2;
  core::CscBuilder pb(2, 2);
  pb.AddEntry(0, 0, 1.0);
  pb.AddEntry(1, 1, 1.0);
  p.p_upper = std::move(pb).Build();
  core::CscBuilder ab(2, 2);
  ab.AddEntry(0, 0, 1.0);
  ab.AddEntry(1, 1, 1.0);
  p.a = std::move(ab).Build();
  p.q = {-1.0, -1.0};
  p.row_lo = {0.0, 0.0};
  p.row_hi = {10.0, 10.0};

  auto s = qp::SolveQpAdmm(p);
  Expect(s.status == core::SolveStatus::kOptimal, "interior QP solves");
  ExpectNear(s.x[0], 1.0, 1e-5, "interior QP x0");
  ExpectNear(s.x[1], 1.0, 1e-5, "interior QP x1");
  ExpectNear(s.objective_value, -1.0, 1e-5, "interior QP objective");
}

// Same quadratic, but the box now cuts the optimum off: x >= 2 forces both
// variables to their lower bound, giving x = (2, 2) and objective
// 1/2(4+4) - 2 - 2 = 0. This is the case that actually exercises the
// projection and the dual update.
void TestActiveBound() {
  core::QpProblem p;
  p.name = "active_bound";
  p.num_cols = 2;
  p.num_rows = 2;
  core::CscBuilder pb(2, 2);
  pb.AddEntry(0, 0, 1.0);
  pb.AddEntry(1, 1, 1.0);
  p.p_upper = std::move(pb).Build();
  core::CscBuilder ab(2, 2);
  ab.AddEntry(0, 0, 1.0);
  ab.AddEntry(1, 1, 1.0);
  p.a = std::move(ab).Build();
  p.q = {-1.0, -1.0};
  p.row_lo = {2.0, 2.0};
  p.row_hi = {10.0, 10.0};

  auto s = qp::SolveQpAdmm(p);
  Expect(s.status == core::SolveStatus::kOptimal, "bound-active QP solves");
  ExpectNear(s.x[0], 2.0, 1e-5, "bound-active QP x0");
  ExpectNear(s.x[1], 2.0, 1e-5, "bound-active QP x1");
  ExpectNear(s.objective_value, 0.0, 1e-5, "bound-active QP objective");
}

// A coupling constraint x0 + x1 = 1 on min 1/2(x0^2 + x1^2). Symmetry puts
// the optimum at (0.5, 0.5) with objective 0.25 -- and the off-diagonal
// coupling is what would expose a mistake in how only P's upper triangle
// is stored and applied.
void TestEqualityCoupling() {
  core::QpProblem p;
  p.name = "coupling";
  p.num_cols = 2;
  p.num_rows = 1;
  core::CscBuilder pb(2, 2);
  pb.AddEntry(0, 0, 1.0);
  pb.AddEntry(1, 1, 1.0);
  p.p_upper = std::move(pb).Build();
  core::CscBuilder ab(1, 2);
  ab.AddEntry(0, 0, 1.0);
  ab.AddEntry(1, 0, 1.0);
  p.a = std::move(ab).Build();
  p.q = {0.0, 0.0};
  p.row_lo = {1.0};
  p.row_hi = {1.0};

  auto s = qp::SolveQpAdmm(p);
  Expect(s.status == core::SolveStatus::kOptimal, "equality-coupled QP solves");
  ExpectNear(s.x[0], 0.5, 1e-5, "equality-coupled QP x0");
  ExpectNear(s.x[1], 0.5, 1e-5, "equality-coupled QP x1");
  ExpectNear(s.objective_value, 0.25, 1e-5, "equality-coupled QP objective");
}

int main() {
  TestUnconstrainedInterior();
  TestActiveBound();
  TestEqualityCoupling();
  if (g_failures) { std::cerr << g_failures << " test(s) failed\n"; return 1; }
  std::cout << "all QP tests passed\n";
  return 0;
}
