#include "models/industrial_mip.hpp"

#include <cmath>
#include <string>
#include <vector>

#include "core/sparse.hpp"

namespace inferno::models {

using core::kInfinity;

namespace {
core::MipProblem Finish(core::MipProblem p,
                         std::vector<std::vector<std::pair<int, double>>>& cols, int rows) {
  p.lp.num_rows = rows;
  core::CscBuilder b(rows, p.lp.num_cols);
  for (int c = 0; c < p.lp.num_cols; ++c) {
    for (const auto& [r, v] : cols[c]) {
      if (v != 0.0) b.AddEntry(c, r, v);
    }
  }
  p.lp.a = std::move(b).Build();
  return p;
}
}  // namespace

core::MipProblem BuildFacilityMip(const FacilityModel& m) {
  const int np = static_cast<int>(m.plant_names.size());
  const int nc = static_cast<int>(m.customer_names.size());
  // Columns: open_p (np binaries), then ship[p][c].
  const int ncols = np + np * nc;
  auto ship = [np, nc](int p, int c) { return np + p * nc + c; };

  core::MipProblem mp;
  mp.lp.name = "FACILITY";
  mp.lp.num_cols = ncols;
  mp.lp.col_lo.assign(ncols, 0.0);
  mp.lp.col_hi.assign(ncols, kInfinity);
  mp.lp.obj.assign(ncols, 0.0);
  mp.lp.col_names.resize(ncols);
  mp.is_integer.assign(ncols, 0);

  for (int p = 0; p < np; ++p) {
    mp.lp.col_hi[p] = 1.0;
    mp.lp.obj[p] = m.fixed_cost[p];
    mp.lp.col_names[p] = "open_" + m.plant_names[p];
    mp.is_integer[p] = 1;
    for (int c = 0; c < nc; ++c) {
      mp.lp.obj[ship(p, c)] = m.transport_cost[p][c];
      mp.lp.col_names[ship(p, c)] = m.plant_names[p] + "->" + m.customer_names[c];
    }
  }

  std::vector<std::vector<std::pair<int, double>>> cols(ncols);
  int row = 0;
  for (int c = 0; c < nc; ++c) {  // demand met
    for (int p = 0; p < np; ++p) cols[ship(p, c)].emplace_back(row, 1.0);
    mp.lp.row_lo.push_back(m.demand[c]);
    mp.lp.row_hi.push_back(kInfinity);
    mp.lp.row_names.push_back("demand_" + m.customer_names[c]);
    ++row;
  }
  for (int p = 0; p < np; ++p) {  // capacity, and the open/ship link
    for (int c = 0; c < nc; ++c) cols[ship(p, c)].emplace_back(row, 1.0);
    cols[p].emplace_back(row, -m.capacity[p]);
    mp.lp.row_lo.push_back(-kInfinity);
    mp.lp.row_hi.push_back(0.0);
    mp.lp.row_names.push_back("cap_" + m.plant_names[p]);
    ++row;
  }
  return Finish(std::move(mp), cols, row);
}

FacilityModel ExampleFacilityModel() {
  FacilityModel m;
  m.plant_names = {"Mangalore", "Chennai", "Kochi"};
  m.fixed_cost = {900.0, 700.0, 800.0};
  m.capacity = {100.0, 80.0, 60.0};
  m.customer_names = {"Bengaluru", "Hyderabad", "Pune"};
  m.demand = {60.0, 50.0, 40.0};
  m.transport_cost = {
      {4.0, 9.0, 12.0},
      {6.0, 5.0, 11.0},
      {5.0, 10.0, 13.0},
  };
  return m;
}

core::MipProblem BuildWorkforceMip(const WorkforceModel& m) {
  const int ns = static_cast<int>(m.shift_names.size());
  const int npat = static_cast<int>(m.pattern_names.size());
  core::MipProblem mp;
  mp.lp.name = "WORKFORCE";
  mp.lp.num_cols = npat;
  mp.lp.col_lo.assign(npat, 0.0);
  mp.lp.col_hi.assign(npat, 50.0);
  mp.lp.obj = m.pattern_cost;
  mp.lp.col_names = m.pattern_names;
  mp.is_integer.assign(npat, 1);  // whole workers only

  std::vector<std::vector<std::pair<int, double>>> cols(npat);
  int row = 0;
  for (int s = 0; s < ns; ++s) {
    for (int k = 0; k < npat; ++k) {
      if (m.covers[s][k] != 0.0) cols[k].emplace_back(row, m.covers[s][k]);
    }
    mp.lp.row_lo.push_back(m.required[s]);
    mp.lp.row_hi.push_back(kInfinity);
    mp.lp.row_names.push_back("cover_" + m.shift_names[s]);
    ++row;
  }
  return Finish(std::move(mp), cols, row);
}

WorkforceModel ExampleWorkforceModel() {
  WorkforceModel m;
  m.shift_names = {"Mon_AM", "Mon_PM", "Tue_AM", "Tue_PM"};
  m.required = {3.0, 5.0, 4.0, 2.0};
  m.shift_cost = {};
  m.max_workers = {};
  // Patterns: each covers two shifts, so shifts are coupled and the
  // relaxation is fractional.
  m.pattern_names = {"P_MonFull", "P_AMs", "P_PMs", "P_TueFull"};
  m.pattern_cost = {100.0, 95.0, 105.0, 98.0};
  m.covers = {
      //  MonFull AMs PMs TueFull
      {1, 1, 0, 0},  // Mon_AM
      {1, 0, 1, 0},  // Mon_PM
      {0, 1, 0, 1},  // Tue_AM
      {0, 0, 1, 1},  // Tue_PM
  };
  return m;
}

core::MipProblem BuildLotSizingMip(const LotSizingModel& m) {
  const int T = static_cast<int>(m.demand.size());
  // Columns: produce_t (T), setup_t (T binary), stock_t (T)
  const int ncols = 3 * T;
  auto prod = [](int t) { return t; };
  auto setup = [T](int t) { return T + t; };
  auto stock = [T](int t) { return 2 * T + t; };

  core::MipProblem mp;
  mp.lp.name = "LOTSIZING";
  mp.lp.num_cols = ncols;
  mp.lp.col_lo.assign(ncols, 0.0);
  mp.lp.col_hi.assign(ncols, kInfinity);
  mp.lp.obj.assign(ncols, 0.0);
  mp.lp.col_names.resize(ncols);
  mp.is_integer.assign(ncols, 0);
  for (int t = 0; t < T; ++t) {
    mp.lp.col_hi[prod(t)] = m.capacity;
    mp.lp.obj[prod(t)] = m.unit_cost;
    mp.lp.col_names[prod(t)] = "produce_t" + std::to_string(t);
    mp.lp.col_hi[setup(t)] = 1.0;
    mp.lp.obj[setup(t)] = m.setup_cost;
    mp.lp.col_names[setup(t)] = "setup_t" + std::to_string(t);
    mp.is_integer[setup(t)] = 1;
    mp.lp.obj[stock(t)] = m.holding_cost;
    mp.lp.col_names[stock(t)] = "stock_t" + std::to_string(t);
  }

  std::vector<std::vector<std::pair<int, double>>> cols(ncols);
  int row = 0;
  // Inventory balance: stock_{t-1} + produce_t - demand_t = stock_t
  for (int t = 0; t < T; ++t) {
    cols[prod(t)].emplace_back(row, 1.0);
    if (t > 0) cols[stock(t - 1)].emplace_back(row, 1.0);
    cols[stock(t)].emplace_back(row, -1.0);
    mp.lp.row_lo.push_back(m.demand[t]);
    mp.lp.row_hi.push_back(m.demand[t]);
    mp.lp.row_names.push_back("balance_t" + std::to_string(t));
    ++row;
  }
  // produce_t <= capacity * setup_t
  for (int t = 0; t < T; ++t) {
    cols[prod(t)].emplace_back(row, 1.0);
    cols[setup(t)].emplace_back(row, -m.capacity);
    mp.lp.row_lo.push_back(-kInfinity);
    mp.lp.row_hi.push_back(0.0);
    mp.lp.row_names.push_back("link_t" + std::to_string(t));
    ++row;
  }
  return Finish(std::move(mp), cols, row);
}

LotSizingModel ExampleLotSizingModel() {
  LotSizingModel m;
  m.demand = {40.0, 20.0, 60.0, 10.0};
  m.unit_cost = 2.0;
  m.setup_cost = 150.0;   // high enough that batching beats producing every period
  m.holding_cost = 1.0;
  m.capacity = 100.0;
  return m;
}

core::MipProblem BuildUnitCommitmentMip(const UnitCommitmentModel& m) {
  const int ng = static_cast<int>(m.unit_names.size());
  const int nt = static_cast<int>(m.load.size());
  // Columns: p[g][t] then on[g][t] (binary)
  const int ncols = 2 * ng * nt;
  auto pcol = [nt](int g, int t) { return g * nt + t; };
  auto oncol = [ng, nt](int g, int t) { return ng * nt + g * nt + t; };

  core::MipProblem mp;
  mp.lp.name = "UNIT_COMMITMENT";
  mp.lp.num_cols = ncols;
  mp.lp.col_lo.assign(ncols, 0.0);
  mp.lp.col_hi.assign(ncols, kInfinity);
  mp.lp.obj.assign(ncols, 0.0);
  mp.lp.col_names.resize(ncols);
  mp.is_integer.assign(ncols, 0);
  for (int g = 0; g < ng; ++g) {
    for (int t = 0; t < nt; ++t) {
      mp.lp.col_hi[pcol(g, t)] = m.pmax[g];
      mp.lp.obj[pcol(g, t)] = m.cost_per_mw[g];
      mp.lp.col_names[pcol(g, t)] = m.unit_names[g] + "_p" + std::to_string(t);
      mp.lp.col_hi[oncol(g, t)] = 1.0;
      mp.lp.obj[oncol(g, t)] = m.startup_cost[g];
      mp.lp.col_names[oncol(g, t)] = m.unit_names[g] + "_on" + std::to_string(t);
      mp.is_integer[oncol(g, t)] = 1;
    }
  }

  std::vector<std::vector<std::pair<int, double>>> cols(ncols);
  int row = 0;
  for (int t = 0; t < nt; ++t) {  // power balance
    for (int g = 0; g < ng; ++g) cols[pcol(g, t)].emplace_back(row, 1.0);
    mp.lp.row_lo.push_back(m.load[t]);
    mp.lp.row_hi.push_back(m.load[t]);
    mp.lp.row_names.push_back("balance_t" + std::to_string(t));
    ++row;
  }
  for (int g = 0; g < ng; ++g) {
    for (int t = 0; t < nt; ++t) {
      // p <= pmax * on
      cols[pcol(g, t)].emplace_back(row, 1.0);
      cols[oncol(g, t)].emplace_back(row, -m.pmax[g]);
      mp.lp.row_lo.push_back(-kInfinity);
      mp.lp.row_hi.push_back(0.0);
      mp.lp.row_names.push_back("ub_" + m.unit_names[g] + "_t" + std::to_string(t));
      ++row;
      // p >= pmin * on  — minimum stable generation, enforced ONLY when
      // committed. This is the constraint the LP relaxation cannot express
      // and the reason unit commitment is a MILP.
      cols[pcol(g, t)].emplace_back(row, 1.0);
      cols[oncol(g, t)].emplace_back(row, -m.pmin[g]);
      mp.lp.row_lo.push_back(0.0);
      mp.lp.row_hi.push_back(kInfinity);
      mp.lp.row_names.push_back("lb_" + m.unit_names[g] + "_t" + std::to_string(t));
      ++row;
    }
  }
  return Finish(std::move(mp), cols, row);
}

UnitCommitmentModel ExampleUnitCommitmentModel() {
  UnitCommitmentModel m;
  m.unit_names = {"Coal", "CCGT", "Peaker"};
  m.cost_per_mw = {21.0, 34.0, 68.0};
  m.startup_cost = {400.0, 250.0, 60.0};
  m.pmin = {200.0, 100.0, 0.0};
  m.pmax = {700.0, 500.0, 400.0};
  m.load = {700.0, 950.0, 1250.0, 820.0};
  return m;
}

}  // namespace inferno::models
