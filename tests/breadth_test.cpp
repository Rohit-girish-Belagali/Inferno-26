// Breadth model tests. As with the refinery model, these recompute the
// structural properties from the raw solution rather than trusting the LP
// — a sign slip in a constraint would otherwise pass unnoticed.
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "checker/checker.hpp"
#include "models/breadth_models.hpp"
#include "simplex/revised_simplex.hpp"

namespace {
int g_failures = 0;
void Expect(bool c, const std::string& m) {
  if (!c) { std::cerr << "FAIL: " << m << "\n"; ++g_failures; }
  else { std::cout << "ok: " << m << "\n"; }
}
}  // namespace

using namespace inferno;

void TestTransportation() {
  models::TransportationModel m = models::ExampleTransportationModel();
  core::LpProblem lp = models::BuildTransportationLp(m);
  core::Solution s = simplex::SolveRevised(lp);
  Expect(s.status == core::SolveStatus::kOptimal, "transportation LP solves");
  if (s.status != core::SolveStatus::kOptimal) return;
  Expect(checker::VerifySolution(lp, s).passed, "transportation passes the independent checker");

  const int ns = static_cast<int>(m.source_names.size());
  const int nd = static_cast<int>(m.sink_names.size());
  for (int i = 0; i < ns; ++i) {
    double ship = 0.0;
    for (int j = 0; j < nd; ++j) ship += s.x[i * nd + j];
    Expect(ship <= m.supply[i] + 1e-6, m.source_names[i] + " ships within its supply");
  }
  for (int j = 0; j < nd; ++j) {
    double recv = 0.0;
    for (int i = 0; i < ns; ++i) recv += s.x[i * nd + j];
    Expect(recv >= m.demand[j] - 1e-6, m.sink_names[j] + " receives its demand");
  }
  for (int k = 0; k < lp.num_cols; ++k) {
    if (s.x[k] < -1e-6) { Expect(false, "no negative shipments"); return; }
  }
  Expect(true, "no negative shipments");

  // Cheapest-source sanity: Bengaluru is cheapest from Mangalore (12 vs 17
  // and 14), so an optimal plan must use that lane. A model with a flipped
  // cost sign would still be feasible, so this checks the objective is
  // actually being minimised rather than merely satisfied.
  Expect(s.x[0 * nd + 0] > 1e-6, "the cheapest lane into Bengaluru is used");
}

void TestEconomicDispatch() {
  models::DispatchModel m = models::ExampleDispatchModel();
  core::LpProblem lp = models::BuildEconomicDispatchLp(m);
  core::Solution s = simplex::SolveRevised(lp);
  Expect(s.status == core::SolveStatus::kOptimal, "economic dispatch LP solves");
  if (s.status != core::SolveStatus::kOptimal) return;
  Expect(checker::VerifySolution(lp, s).passed, "dispatch passes the independent checker");

  const int ng = static_cast<int>(m.generators.size());
  const int nt = static_cast<int>(m.load.size());
  for (int t = 0; t < nt; ++t) {
    double gen = 0.0;
    for (int g = 0; g < ng; ++g) gen += s.x[g * nt + t];
    Expect(std::abs(gen - m.load[t]) < 1e-6,
           "period " + std::to_string(t) + " generation meets load exactly");
  }
  for (int g = 0; g < ng; ++g) {
    for (int t = 0; t < nt; ++t) {
      double p = s.x[g * nt + t];
      Expect(p >= m.generators[g].pmin - 1e-6 && p <= m.generators[g].pmax + 1e-6,
             m.generators[g].name + " within limits at t" + std::to_string(t));
    }
    for (int t = 0; t + 1 < nt; ++t) {
      double d = std::abs(s.x[g * nt + t + 1] - s.x[g * nt + t]);
      Expect(d <= m.generators[g].ramp + 1e-6,
             m.generators[g].name + " respects its ramp limit at t" + std::to_string(t));
    }
  }

  // Merit order: the peaker is the most expensive unit, so it should sit
  // at zero in the lightest period. If it runs when it need not, the
  // objective is not being minimised.
  Expect(s.x[2 * nt + 0] < 1e-6, "the expensive peaker is off in the lightest period");
}

int main() {
  TestTransportation();
  TestEconomicDispatch();
  if (g_failures) { std::cerr << g_failures << " test(s) failed\n"; return 1; }
  std::cout << "all breadth model tests passed\n";
  return 0;
}
