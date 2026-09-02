#pragma once

#include <string>

#include "core/lp_problem.hpp"

namespace inferno::io {

// Writes a solved core::Solution alongside the problem it was solved for to
// a simple, greppable text format: one line per column, then one line per
// row, then the objective value and status.
void WriteSolution(const std::string& path, const core::LpProblem& problem,
                    const core::Solution& solution);

std::string StatusToString(core::SolveStatus status);

}  // namespace inferno::io
