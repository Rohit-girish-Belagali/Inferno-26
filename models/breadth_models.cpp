#include "models/breadth_models.hpp"

#include <cmath>
#include <vector>

#include "core/sparse.hpp"

namespace inferno::models {

using core::kInfinity;

core::LpProblem BuildTransportationLp(const TransportationModel& m) {
  const int ns = static_cast<int>(m.source_names.size());
  const int nd = static_cast<int>(m.sink_names.size());
  const int ncols = ns * nd;
  auto col = [nd](int i, int j) { return i * nd + j; };

  core::LpProblem lp;
  lp.name = "TRANSPORT";
  lp.num_cols = ncols;
  lp.col_lo.assign(ncols, 0.0);
  lp.col_hi.assign(ncols, kInfinity);
  lp.obj.assign(ncols, 0.0);
  lp.col_names.resize(ncols);
  for (int i = 0; i < ns; ++i) {
    for (int j = 0; j < nd; ++j) {
      lp.obj[col(i, j)] = m.cost[i][j];
      lp.col_names[col(i, j)] = m.source_names[i] + "->" + m.sink_names[j];
    }
  }

  std::vector<std::vector<std::pair<int, double>>> cols(ncols);
  int row = 0;
  for (int i = 0; i < ns; ++i) {  // supply
    for (int j = 0; j < nd; ++j) cols[col(i, j)].emplace_back(row, 1.0);
    lp.row_lo.push_back(-kInfinity);
    lp.row_hi.push_back(m.supply[i]);
    lp.row_names.push_back("SUP_" + m.source_names[i]);
    ++row;
  }
  for (int j = 0; j < nd; ++j) {  // demand
    for (int i = 0; i < ns; ++i) cols[col(i, j)].emplace_back(row, 1.0);
    lp.row_lo.push_back(m.demand[j]);
    lp.row_hi.push_back(kInfinity);
    lp.row_names.push_back("DEM_" + m.sink_names[j]);
    ++row;
  }

  lp.num_rows = row;
  core::CscBuilder b(row, ncols);
  for (int c = 0; c < ncols; ++c) {
    for (const auto& [r, v] : cols[c]) b.AddEntry(c, r, v);
  }
  lp.a = std::move(b).Build();
  return lp;
}

TransportationModel ExampleTransportationModel() {
  TransportationModel m;
  m.source_names = {"Mangalore", "Kochi", "Chennai"};
  m.sink_names = {"Bengaluru", "Hyderabad", "Pune", "Nagpur"};
  m.supply = {1800.0, 1200.0, 1500.0};
  m.demand = {1400.0, 1000.0, 1100.0, 700.0};  // 4200 total against 4500 supply
  // Roughly distance-proportional freight cost per unit.
  m.cost = {
      {12.0, 22.0, 31.0, 38.0},
      {17.0, 26.0, 34.0, 42.0},
      {14.0, 15.0, 29.0, 30.0},
  };
  return m;
}

core::LpProblem BuildEconomicDispatchLp(const DispatchModel& m) {
  const int ng = static_cast<int>(m.generators.size());
  const int nt = static_cast<int>(m.load.size());
  const int ncols = ng * nt;
  auto col = [nt](int g, int t) { return g * nt + t; };

  core::LpProblem lp;
  lp.name = "ECON_DISPATCH";
  lp.num_cols = ncols;
  lp.col_lo.resize(ncols);
  lp.col_hi.resize(ncols);
  lp.obj.resize(ncols);
  lp.col_names.resize(ncols);
  for (int g = 0; g < ng; ++g) {
    for (int t = 0; t < nt; ++t) {
      lp.col_lo[col(g, t)] = m.generators[g].pmin;
      lp.col_hi[col(g, t)] = m.generators[g].pmax;
      lp.obj[col(g, t)] = m.generators[g].cost_per_mw;
      lp.col_names[col(g, t)] = m.generators[g].name + "_t" + std::to_string(t);
    }
  }

  std::vector<std::vector<std::pair<int, double>>> cols(ncols);
  int row = 0;
  for (int t = 0; t < nt; ++t) {  // power balance, an equality
    for (int g = 0; g < ng; ++g) cols[col(g, t)].emplace_back(row, 1.0);
    lp.row_lo.push_back(m.load[t]);
    lp.row_hi.push_back(m.load[t]);
    lp.row_names.push_back("BAL_t" + std::to_string(t));
    ++row;
  }
  // Ramp limits, as a two-sided range on the difference between periods:
  // -ramp <= p_{t+1} - p_t <= ramp. One row rather than two, since this
  // representation carries both bounds natively.
  for (int g = 0; g < ng; ++g) {
    for (int t = 0; t + 1 < nt; ++t) {
      cols[col(g, t + 1)].emplace_back(row, 1.0);
      cols[col(g, t)].emplace_back(row, -1.0);
      lp.row_lo.push_back(-m.generators[g].ramp);
      lp.row_hi.push_back(m.generators[g].ramp);
      lp.row_names.push_back("RAMP_" + m.generators[g].name + "_t" + std::to_string(t));
      ++row;
    }
  }

  lp.num_rows = row;
  core::CscBuilder b(row, ncols);
  for (int c = 0; c < ncols; ++c) {
    for (const auto& [r, v] : cols[c]) b.AddEntry(c, r, v);
  }
  lp.a = std::move(b).Build();
  return lp;
}

DispatchModel ExampleDispatchModel() {
  DispatchModel m;
  // A cheap inflexible baseload unit, a mid-merit unit, and an expensive
  // fast peaker -- the classic merit order, and the ramp limits are what
  // force the expensive unit on during the evening peak even though the
  // cheap one has headroom left.
  m.generators = {
      {"Coal_Base", 2100.0, 200.0, 700.0, 80.0},
      {"Gas_CCGT", 3400.0, 100.0, 500.0, 200.0},
      {"Gas_Peaker", 6800.0, 0.0, 400.0, 400.0},
  };
  m.load = {700.0, 780.0, 950.0, 1250.0, 1100.0, 820.0};
  return m;
}

}  // namespace inferno::models
