#pragma once

#include <string>
#include <vector>

#include "core/lp_problem.hpp"

namespace inferno::models {

// The industrial and stress-test problem battery for PS 26119.
//
// The problem statement names specific sectors — refinery scheduling,
// crude blending, production planning, logistics, power dispatch,
// transportation, supply chain — and specific stress conditions:
// degeneracy, ill-conditioning, and scale into the millions of variables.
// These are the ones expressible as pure LP. Anything requiring integer
// decisions is NOT here, because this project has no MILP solver; see
// STATUS.md, which states that gap rather than working around it.

// --- Level 1: production planning -------------------------------------
// Maximise profit over a product mix subject to shared resource limits.
// The canonical first LP, and the one every reviewer can verify by hand.
struct ProductionModel {
  std::vector<std::string> product_names;
  std::vector<double> profit;                    // per unit
  std::vector<std::string> resource_names;
  std::vector<double> available;                 // per resource
  std::vector<std::vector<double>> usage;        // usage[resource][product]
};
core::LpProblem BuildProductionLp(const ProductionModel& m);
ProductionModel ExampleProductionModel();

// --- Level 1: diet ----------------------------------------------------
// Cheapest food basket meeting nutritional minima — Stigler's problem,
// and historically one of the first LPs ever solved.
struct DietModel {
  std::vector<std::string> food_names;
  std::vector<double> cost;                      // per unit of food
  std::vector<std::string> nutrient_names;
  std::vector<double> min_required;              // per nutrient
  std::vector<double> max_allowed;               // per nutrient, +inf if none
  std::vector<std::vector<double>> content;      // content[nutrient][food]
};
core::LpProblem BuildDietLp(const DietModel& m);
DietModel ExampleDietModel();

// --- Stress: degeneracy -----------------------------------------------
// A transportation-style LP built so that a large number of distinct
// bases give the identical objective value. Degeneracy is what makes a
// naive simplex stall: the ratio test returns a zero-length step and the
// method can pivot indefinitely without improving. `blocks` controls the
// size; the variable count grows roughly as blocks^2.
core::LpProblem BuildDegenerateLp(int blocks);

// --- Stress: ill-conditioning -----------------------------------------
// Coefficients spanning `decades` orders of magnitude, from 10^-k to
// 10^+k. This is where an unscaled solver loses digits in elimination and
// starts reporting confident nonsense; the PS calls for exactly this test.
core::LpProblem BuildIllConditionedLp(int n, int decades);

// --- Stress: scale ----------------------------------------------------
// A sparse LP with `cols` variables and `rows` constraints, roughly
// `nnz_per_col` nonzeros per column, generated deterministically from
// `seed` so a run is reproducible. Used to characterise how the engine
// behaves as the model grows toward the millions of variables the PS
// names.
core::LpProblem BuildLargeSparseLp(int rows, int cols, int nnz_per_col, unsigned seed);

}  // namespace inferno::models
