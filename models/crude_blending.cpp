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


core::MipProblem BuildCrudeBlendingMip(const RefineryMipModel& m) {
  // Start from the LP, then append one activation binary per crude and a
  // linking row per crude. Reusing the LP builder verbatim is deliberate:
  // it guarantees the blending and quality rows are byte-identical to the
  // continuous model, so the MILP is provably the same problem plus the
  // discrete decision, not a re-derivation that might differ subtly.
  core::LpProblem lp = BuildCrudeBlendingLp(m.blending);
  const int nc = static_cast<int>(m.blending.crudes.size());
  const int np = static_cast<int>(m.blending.products.size());
  const int base_cols = lp.num_cols;
  const int base_rows = lp.num_rows;

  core::MipProblem mp;
  mp.lp = lp;
  mp.lp.name = "CRUDE_BLEND_MIP";
  mp.lp.num_cols = base_cols + nc;
  mp.lp.col_lo.resize(mp.lp.num_cols, 0.0);
  mp.lp.col_hi.resize(mp.lp.num_cols, 1.0);
  mp.lp.obj.resize(mp.lp.num_cols, 0.0);
  mp.lp.col_names.resize(mp.lp.num_cols);
  mp.is_integer.assign(mp.lp.num_cols, 0);
  for (int i = 0; i < nc; ++i) {
    mp.lp.obj[base_cols + i] = m.activation_cost[i];
    mp.lp.col_names[base_cols + i] = "use_" + m.blending.crudes[i].name;
    mp.is_integer[base_cols + i] = 1;
  }

  // Rebuild the matrix with the original entries plus the link rows:
  //   sum_j x_ij - availability_i * use_i <= 0
  // so a crude can only supply volume when its binary is set.
  std::vector<std::vector<std::pair<int, double>>> cols(mp.lp.num_cols);
  for (int c = 0; c < base_cols; ++c) {
    for (int p = lp.a.col_ptr[c]; p < lp.a.col_ptr[c + 1]; ++p) {
      cols[c].emplace_back(lp.a.row_idx[p], lp.a.values[p]);
    }
  }
  int row = base_rows;
  if (m.throughput_capacity > 0.0) {
    for (int i = 0; i < nc; ++i) {
      for (int j = 0; j < np; ++j) cols[i * np + j].emplace_back(row, 1.0);
    }
    mp.lp.row_lo.push_back(-kInfinity);
    mp.lp.row_hi.push_back(m.throughput_capacity);
    mp.lp.row_names.push_back("throughput");
    ++row;
  }
  for (int i = 0; i < nc; ++i) {
    for (int j = 0; j < np; ++j) cols[i * np + j].emplace_back(row, 1.0);
    cols[base_cols + i].emplace_back(row, -m.blending.crudes[i].availability);
    mp.lp.row_lo.push_back(-kInfinity);
    mp.lp.row_hi.push_back(0.0);
    mp.lp.row_names.push_back("use_link_" + m.blending.crudes[i].name);
    ++row;
  }
  mp.lp.num_rows = row;
  core::CscBuilder b(row, mp.lp.num_cols);
  for (int c = 0; c < mp.lp.num_cols; ++c) {
    for (const auto& [r, v] : cols[c]) {
      if (v != 0.0) b.AddEntry(c, r, v);
    }
  }
  mp.lp.a = std::move(b).Build();
  return mp;
}

RefineryMipModel ExampleRefineryMipModel() {
  RefineryMipModel m;
  m.blending = ExampleRefineryModel();
  // Commitment costs sized so activating every crude is NOT automatically
  // optimal. This mattered: the first values used here were low enough
  // that all four binaries came out at 1 in the relaxation, the model
  // solved at zero nodes, and the discrete decision was vacuous -- the
  // test passed while exercising nothing. Condensate contributes 12,000
  // bbl at a marginal value near 33.8/bbl, so a commitment cost above
  // roughly 405,000 makes buying it a genuinely open question and forces
  // the search to decide rather than to agree with the relaxation.
  m.activation_cost = {120000.0, 90000.0, 60000.0, 600000.0};
  // Below the 137,000 bbl of crude available, so the refinery must choose
  // which barrels to buy rather than taking all of them.
  m.throughput_capacity = 110000.0;
  return m;
}

}  // namespace inferno::models
