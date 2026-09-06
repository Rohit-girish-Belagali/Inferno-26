#pragma once

#include <string>
#include <vector>

#include "core/lp_problem.hpp"

namespace inferno::models {

// Crude blending — the refinery model BUILD_PLAN_V2.md Phase 5 asks for
// first, and the reason it does: "Generic benchmark numbers will not move
// MRPL. A working crude blending model will."
//
// This is the classical linear blending formulation. A refinery buys
// several crudes, each with its own cost, availability and quality
// properties (sulfur, density, and so on), and allocates them across
// finished products, each of which has a demand to meet, a selling price,
// and specification limits on those same properties.
//
// Decision variables x[i][j] — volume of crude i blended into product j.
//
//   maximise   sum_j price_j * (sum_i x_ij)  -  sum_i cost_i * (sum_j x_ij)
//   subject to sum_j x_ij <= availability_i                     (supply)
//              sum_i x_ij >= demand_j                           (demand)
//              quality specs, per product and property
//
// The quality constraint is the part worth stating carefully, because the
// obvious way to write it is nonlinear. A property of a blend is its
// volume-weighted average,
//
//     (sum_i prop_ik x_ij) / (sum_i x_ij)  <=  max_jk
//
// which is a ratio of decision variables. Multiplying through by the
// (nonnegative) denominator clears it to
//
//     sum_i (prop_ik - max_jk) * x_ij  <=  0
//
// which is linear, and exactly equivalent whenever the blend is nonempty.
// That transformation is what keeps this an LP rather than requiring the
// MILP machinery this project deliberately did not build (the plan's own
// day-25 kill checkpoint says to narrow scope to LP rather than ship a
// shaky MILP). It holds for properties that blend linearly by volume;
// properties that do not (viscosity, octane, pour point) blend by an
// index rather than directly, so the caller is expected to supply
// already-indexed values — noted here rather than silently assumed.
struct Crude {
  std::string name;
  double cost = 0.0;          // currency per barrel
  double availability = 0.0;  // barrels available this period
  std::vector<double> properties;  // one entry per property, same order everywhere
};

struct Product {
  std::string name;
  double price = 0.0;   // currency per barrel
  double demand = 0.0;  // barrels that must be produced
  // Per-property specification window. Use -infinity / +infinity for an
  // unconstrained side.
  std::vector<double> min_spec;
  std::vector<double> max_spec;
};

struct BlendingModel {
  std::vector<std::string> property_names;
  std::vector<Crude> crudes;
  std::vector<Product> products;
};

// Builds the LP. Column ordering is x[i][j] at index i * products + j, so
// a solution vector can be read back without a lookup table.
core::LpProblem BuildCrudeBlendingLp(const BlendingModel& model);

// A small but realistic instance for demos and tests: four crudes with
// genuinely different sulfur/density profiles and costs, three products
// with real specification windows.
BlendingModel ExampleRefineryModel();

}  // namespace inferno::models
