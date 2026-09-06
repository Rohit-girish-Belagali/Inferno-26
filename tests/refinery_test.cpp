// Phase 5 refinery model test. Checks the crude-blending LP is built
// correctly, solves, and — the part that matters — that the answer is
// physically and economically coherent, not merely that the solver
// returned something.
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "checker/checker.hpp"
#include "models/crude_blending.hpp"
#include "simplex/revised_simplex.hpp"

namespace {
int g_failures = 0;
void Expect(bool cond, const std::string& msg) {
  if (!cond) { std::cerr << "FAIL: " << msg << "\n"; ++g_failures; }
  else { std::cout << "ok: " << msg << "\n"; }
}
}  // namespace

int main() {
  using namespace inferno;
  models::BlendingModel model = models::ExampleRefineryModel();
  core::LpProblem lp = models::BuildCrudeBlendingLp(model);

  const int nc = static_cast<int>(model.crudes.size());
  const int np = static_cast<int>(model.products.size());
  Expect(lp.num_cols == nc * np, "one column per (crude, product) pair");

  core::Solution sol = simplex::SolveRevised(lp);
  Expect(sol.status == core::SolveStatus::kOptimal, "refinery blending LP solves to optimality");
  if (sol.status != core::SolveStatus::kOptimal) return 1;

  auto chk = checker::VerifySolution(lp, sol);
  Expect(chk.passed, "refinery solution passes the independent checker: " + chk.message);

  // Supply must never be exceeded.
  for (int i = 0; i < nc; ++i) {
    double used = 0.0;
    for (int j = 0; j < np; ++j) used += sol.x[i * np + j];
    Expect(used <= model.crudes[i].availability + 1e-6,
           model.crudes[i].name + " respects its availability");
  }

  // Demand must be met, and every blend must honour its own spec window —
  // recomputed here from the raw volumes rather than trusted from the LP,
  // since a sign slip in the quality rows would otherwise pass unnoticed.
  for (int j = 0; j < np; ++j) {
    double vol = 0.0;
    std::vector<double> prop(model.property_names.size(), 0.0);
    for (int i = 0; i < nc; ++i) {
      double v = sol.x[i * np + j];
      vol += v;
      for (size_t k = 0; k < prop.size(); ++k) prop[k] += v * model.crudes[i].properties[k];
    }
    Expect(vol >= model.products[j].demand - 1e-6, model.products[j].name + " meets demand");
    for (size_t k = 0; k < prop.size(); ++k) {
      if (vol <= 0.0) continue;
      double avg = prop[k] / vol;
      double lo = model.products[j].min_spec[k], hi = model.products[j].max_spec[k];
      if (std::isfinite(hi)) {
        Expect(avg <= hi + 1e-6,
               model.products[j].name + " " + model.property_names[k] + " within max spec");
      }
      if (std::isfinite(lo)) {
        Expect(avg >= lo - 1e-6,
               model.products[j].name + " " + model.property_names[k] + " within min spec");
      }
    }
  }

  Expect(-sol.objective_value > 0.0, "the optimal plan is profitable");

  if (g_failures) { std::cerr << g_failures << " test(s) failed\n"; return 1; }
  std::cout << "all refinery tests passed\n";
  return 0;
}
