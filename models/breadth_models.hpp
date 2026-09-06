#pragma once

#include <string>
#include <vector>

#include "core/lp_problem.hpp"

namespace inferno::models {

// Breadth models (BUILD_PLAN_V2.md Phase 5). The refinery model shows the
// engine on the problem statement's own domain; these show it is a general
// LP engine and not something shaped around a single use case. Both are
// textbook formulations from open literature, which is the point — the
// data is defensible and a reviewer can check the model against a
// standard reference rather than taking our word for the structure.

// ---------------------------------------------------------------------
// Transportation / supply chain.
//
//     min  sum_ij cost_ij x_ij
//     s.t. sum_j x_ij <= supply_i      for each source
//          sum_i x_ij >= demand_j      for each sink
//          x_ij >= 0
//
// The classical Hitchcock transportation problem. Feasible exactly when
// total supply >= total demand, which BuildTransportationLp does NOT
// silently repair: an infeasible instance should be reported infeasible by
// the solver, not quietly patched by the model builder.
struct TransportationModel {
  std::vector<std::string> source_names;
  std::vector<std::string> sink_names;
  std::vector<double> supply;              // per source
  std::vector<double> demand;              // per sink
  std::vector<std::vector<double>> cost;   // cost[source][sink]
};

// Column j of source i is at index i * sinks + j.
core::LpProblem BuildTransportationLp(const TransportationModel& m);
TransportationModel ExampleTransportationModel();

// ---------------------------------------------------------------------
// Power system economic dispatch.
//
//     min  sum_gt (a_g * p_gt)
//     s.t. sum_g p_gt = load_t             for each period   (power balance)
//          pmin_g <= p_gt <= pmax_g                          (generator limits)
//          |p_g,t+1 - p_gt| <= ramp_g                        (ramp limits)
//
// NAMED HONESTLY. The plan's checklist says "unit commitment", and this is
// NOT that: unit commitment decides which generators are switched on,
// which is a binary decision and therefore a MILP — deliberately out of
// scope here (see the day-25 kill checkpoint). This is the economic
// dispatch problem, the continuous sub-problem that unit commitment solves
// repeatedly once the on/off pattern is fixed. It is a real and widely
// used model in its own right, and calling it unit commitment would be
// claiming a capability this project does not have.
//
// Generators are assumed already committed, so pmin is enforced as a hard
// lower bound. Ramp limits couple consecutive periods, which is what makes
// this more than a sequence of independent single-period problems.
struct Generator {
  std::string name;
  double cost_per_mw = 0.0;
  double pmin = 0.0;
  double pmax = 0.0;
  double ramp = 0.0;  // max change in output between consecutive periods
};

struct DispatchModel {
  std::vector<Generator> generators;
  std::vector<double> load;  // one entry per period
};

// Column for generator g in period t is at index g * periods + t.
core::LpProblem BuildEconomicDispatchLp(const DispatchModel& m);
DispatchModel ExampleDispatchModel();

}  // namespace inferno::models
