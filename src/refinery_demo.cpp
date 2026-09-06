// Refinery crude-blending demo (BUILD_PLAN_V2.md Phase 5). Builds the LP,
// solves it with this project's own simplex, verifies the answer with the
// independent checker, and prints the plan a refinery planner would
// actually read — including the shadow prices, which are usually the part
// they care about most.
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

#include "checker/checker.hpp"
#include "models/crude_blending.hpp"
#include "simplex/revised_simplex.hpp"

using namespace inferno;

int main() {
  models::BlendingModel model = models::ExampleRefineryModel();
  core::LpProblem lp = models::BuildCrudeBlendingLp(model);

  const int nc = static_cast<int>(model.crudes.size());
  const int np = static_cast<int>(model.products.size());

  printf("Refinery crude blending — %d crudes, %d products, %d properties\n",
         nc, np, static_cast<int>(model.property_names.size()));
  printf("LP size: %d rows x %d columns, %d nonzeros\n\n", lp.num_rows, lp.num_cols, lp.a.nnz());

  core::Solution sol = simplex::SolveRevised(lp);
  if (sol.status != core::SolveStatus::kOptimal) {
    printf("solver did not reach optimality (status %d)\n", static_cast<int>(sol.status));
    return 1;
  }

  checker::CheckResult chk = checker::VerifySolution(lp, sol);
  printf("Optimal profit: %.2f per period\n", -sol.objective_value);
  printf("Independent checker: %s\n\n", chk.message.c_str());

  printf("%-12s", "BLEND (bbl)");
  for (int j = 0; j < np; ++j) printf("%12s", model.products[j].name.c_str());
  printf("%12s\n", "used");
  for (int i = 0; i < nc; ++i) {
    printf("%-12s", model.crudes[i].name.c_str());
    double used = 0.0;
    for (int j = 0; j < np; ++j) {
      double v = sol.x[i * np + j];
      used += v;
      printf("%12.0f", v);
    }
    printf("%12.0f  of %.0f available\n", used, model.crudes[i].availability);
  }

  printf("\n%-12s%12s%12s%12s\n", "PRODUCT", "volume", "sulfur%", "API");
  for (int j = 0; j < np; ++j) {
    double vol = 0.0;
    std::vector<double> prop(model.property_names.size(), 0.0);
    for (int i = 0; i < nc; ++i) {
      double v = sol.x[i * np + j];
      vol += v;
      for (size_t k = 0; k < prop.size(); ++k) prop[k] += v * model.crudes[i].properties[k];
    }
    printf("%-12s%12.0f", model.products[j].name.c_str(), vol);
    for (size_t k = 0; k < prop.size(); ++k) printf("%12.3f", vol > 0 ? prop[k] / vol : 0.0);
    printf("   (demand %.0f)\n", model.products[j].demand);
  }

  // Shadow prices: what one more barrel of a binding resource is worth.
  // This is the output a planner actually uses to decide what to buy next,
  // so it is worth surfacing rather than leaving buried in the dual vector.
  //
  // Sign: the LP minimises negated profit, so its duals carry the opposite
  // sign to the quantity a planner wants. Reported here in PROFIT terms —
  // positive means one more unit of that constraint's resource is worth
  // that much, negative means the constraint is costing that much to
  // satisfy.
  printf("\nMARGINAL VALUES, in profit terms (binding constraints only)\n");
  bool any = false;
  for (int i = 0; i < lp.num_rows; ++i) {
    if (std::abs(sol.y[i]) < 1e-7) continue;
    printf("  %-28s %+10.3f per unit\n", lp.row_names[i].c_str(), -sol.y[i]);
    any = true;
  }
  if (!any) printf("  (none binding)\n");
  return chk.passed ? 0 : 2;
}
