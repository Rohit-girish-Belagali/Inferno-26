#pragma once

#include <string>
#include <vector>

#include "core/mip_problem.hpp"

namespace inferno::models {

// Industrial MILP models — the problem classes PS 26119 names that need
// integer decisions. Each is small enough that its structure can be read
// and checked by hand, and each is built so the LP relaxation is genuinely
// fractional: a model whose relaxation is already integral would exercise
// no branching and prove nothing about the search.

// --- Plant / facility selection --------------------------------------
// Open a subset of plants (binary) and ship from them to customers
// (continuous), paying fixed opening costs plus per-unit transport, and
// meeting every customer's demand without exceeding any plant's capacity.
// This is the canonical fixed-charge structure: the binary is linked to
// the continuous flow by a capacity constraint, which is what makes the
// relaxation fractional -- it opens a plant "partially".
struct FacilityModel {
  std::vector<std::string> plant_names;
  std::vector<double> fixed_cost;
  std::vector<double> capacity;
  std::vector<std::string> customer_names;
  std::vector<double> demand;
  std::vector<std::vector<double>> transport_cost;  // [plant][customer]
};
core::MipProblem BuildFacilityMip(const FacilityModel& m);
FacilityModel ExampleFacilityModel();

// --- Workforce shift scheduling ---------------------------------------
// Assign whole workers (integer) to shifts so every shift meets its
// minimum staffing at least cost. Workers are indivisible, which is the
// entire reason this is not an LP.
struct WorkforceModel {
  std::vector<std::string> shift_names;
  std::vector<double> required;      // minimum workers per shift
  std::vector<double> shift_cost;    // cost per worker assigned to that shift
  std::vector<double> max_workers;   // availability cap per shift
  // Each pattern covers several shifts at once (a worker on a pattern
  // covers every shift it includes), which is what couples the shifts.
  std::vector<std::string> pattern_names;
  std::vector<double> pattern_cost;
  std::vector<std::vector<double>> covers;  // [shift][pattern] 1 if covered
};
core::MipProblem BuildWorkforceMip(const WorkforceModel& m);
WorkforceModel ExampleWorkforceModel();

// --- Multi-period production planning with setup costs ----------------
// Produce over several periods against demand, holding inventory between
// them, paying a fixed setup cost in any period with positive production.
// The setup binary linked to production by a big-M capacity row is the
// classic lot-sizing structure.
struct LotSizingModel {
  std::vector<double> demand;        // per period
  double unit_cost = 0.0;
  double setup_cost = 0.0;
  double holding_cost = 0.0;
  double capacity = 0.0;             // max production per period
};
core::MipProblem BuildLotSizingMip(const LotSizingModel& m);
LotSizingModel ExampleLotSizingModel();

// --- Power unit commitment --------------------------------------------
// The MILP the economic-dispatch LP could not express: which generators
// are ON in each period (binary), with minimum stable generation enforced
// only when a unit is committed. That on/off decision is precisely what
// made the LP version a relaxation rather than the real problem.
struct UnitCommitmentModel {
  std::vector<std::string> unit_names;
  std::vector<double> cost_per_mw;
  std::vector<double> startup_cost;  // charged per period a unit is on
  std::vector<double> pmin;
  std::vector<double> pmax;
  std::vector<double> load;          // per period
};
core::MipProblem BuildUnitCommitmentMip(const UnitCommitmentModel& m);
UnitCommitmentModel ExampleUnitCommitmentModel();

}  // namespace inferno::models
