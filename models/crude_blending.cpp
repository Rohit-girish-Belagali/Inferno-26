#include "models/crude_blending.hpp"

#include <cmath>

#include "core/sparse.hpp"

namespace inferno::models {

using core::kInfinity;

core::LpProblem BuildCrudeBlendingLp(const BlendingModel& model) {
  const int nc = static_cast<int>(model.crudes.size());
  const int np = static_cast<int>(model.products.size());
  const int nk = static_cast<int>(model.property_names.size());
  const int ncols = nc * np;

  // Rows: one supply row per crude, one demand row per product, and up to
  // two spec rows per (product, property) pair.
  core::LpProblem lp;
  lp.name = "CRUDE_BLEND";
  lp.num_cols = ncols;
  lp.col_lo.assign(ncols, 0.0);
  lp.col_hi.assign(ncols, kInfinity);
  lp.obj.assign(ncols, 0.0);
  lp.col_names.resize(ncols);

  auto col = [np](int i, int j) { return i * np + j; };

  // This project minimises, so the profit objective is negated: each
  // barrel of crude i sold as product j earns price_j and costs cost_i.
  for (int i = 0; i < nc; ++i) {
    for (int j = 0; j < np; ++j) {
      lp.obj[col(i, j)] = -(model.products[j].price - model.crudes[i].cost);
      lp.col_names[col(i, j)] = model.crudes[i].name + "->" + model.products[j].name;
    }
  }

  core::CscBuilder builder(0, ncols);  // row count patched in below
  std::vector<std::vector<std::pair<int, double>>> cols(ncols);
  int row = 0;
  auto add = [&](int c, double v) { cols[c].emplace_back(row, v); };

  // --- Supply: sum_j x_ij <= availability_i ---
  for (int i = 0; i < nc; ++i) {
    for (int j = 0; j < np; ++j) add(col(i, j), 1.0);
    lp.row_lo.push_back(-kInfinity);
    lp.row_hi.push_back(model.crudes[i].availability);
    lp.row_names.push_back("SUP_" + model.crudes[i].name);
    ++row;
  }

  // --- Demand: sum_i x_ij >= demand_j ---
  for (int j = 0; j < np; ++j) {
    for (int i = 0; i < nc; ++i) add(col(i, j), 1.0);
    lp.row_lo.push_back(model.products[j].demand);
    lp.row_hi.push_back(kInfinity);
    lp.row_names.push_back("DEM_" + model.products[j].name);
    ++row;
  }

  // --- Quality: sum_i (prop_ik - spec_jk) x_ij  <= 0  (upper spec)
  //              sum_i (prop_ik - spec_jk) x_ij  >= 0  (lower spec) ---
  for (int j = 0; j < np; ++j) {
    for (int k = 0; k < nk; ++k) {
      double hi = model.products[j].max_spec[k];
      if (std::isfinite(hi)) {
        for (int i = 0; i < nc; ++i) add(col(i, j), model.crudes[i].properties[k] - hi);
        lp.row_lo.push_back(-kInfinity);
        lp.row_hi.push_back(0.0);
        lp.row_names.push_back("MAX_" + model.products[j].name + "_" + model.property_names[k]);
        ++row;
      }
      double lo = model.products[j].min_spec[k];
      if (std::isfinite(lo)) {
        for (int i = 0; i < nc; ++i) add(col(i, j), model.crudes[i].properties[k] - lo);
        lp.row_lo.push_back(0.0);
        lp.row_hi.push_back(kInfinity);
        lp.row_names.push_back("MIN_" + model.products[j].name + "_" + model.property_names[k]);
        ++row;
      }
    }
  }

  lp.num_rows = row;
  core::CscBuilder b(row, ncols);
  for (int c = 0; c < ncols; ++c) {
    for (const auto& [r, v] : cols[c]) {
      if (v != 0.0) b.AddEntry(c, r, v);
    }
  }
  lp.a = std::move(b).Build();
  return lp;
}

BlendingModel ExampleRefineryModel() {
  BlendingModel m;
  // Sulfur in weight %, density as API gravity. Both blend linearly by
  // volume to a good approximation, which is what keeps this an LP.
  m.property_names = {"sulfur_pct", "api_gravity"};

  // Representative of the sweet/sour spread a coastal Indian refinery
  // actually chooses between: light sweet crude is expensive and clean,
  // heavy sour is cheap and needs blending down to meet sulfur specs.
  m.crudes = {
      {"LightSweet", 78.0, 30000.0, {0.20, 38.0}},
      {"MediumSour", 68.0, 45000.0, {1.60, 31.0}},
      {"HeavySour", 61.0, 50000.0, {2.90, 22.0}},
      {"Condensate", 82.0, 12000.0, {0.05, 52.0}},
  };

  // These are specs on the CRUDE SLATE fed to each processing mode, not on
  // finished product — a distinction worth being explicit about, because
  // finished-diesel sulfur limits are orders of magnitude tighter and are
  // met by downstream hydrotreating, not by blending crude. Written at
  // finished-product levels this model is simply infeasible: capping the
  // diesel slate at 0.50% sulfur admits at most ~21,600 bbl against a
  // 40,000 bbl demand, which is how the first version of this data was
  // caught.
  m.products = {
      // Diesel-mode slate: moderate sulfur ceiling, needs a minimum gravity.
      {"Diesel", 96.0, 40000.0, {-kInfinity, 30.0}, {1.20, kInfinity}},
      // Gasoline-mode slate: lighter still, tighter sulfur allowance.
      {"Gasoline", 104.0, 25000.0, {-kInfinity, 34.0}, {0.80, kInfinity}},
      // Fuel-oil mode: the sink for heavy sour barrels, loose specs.
      {"FuelOil", 71.0, 30000.0, {-kInfinity, -kInfinity}, {3.00, kInfinity}},
  };
  return m;
}

}  // namespace inferno::models
