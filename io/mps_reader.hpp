#pragma once

#include <string>

#include "core/lp_problem.hpp"

namespace inferno::io {

// Reads an MPS file (free-form: whitespace-delimited fields, which also
// correctly parses the vast majority of fixed-form files including the
// entire Netlib set) into a core::LpProblem. Throws std::runtime_error on
// malformed input.
//
// Supported sections: NAME, ROWS, COLUMNS (with MARKER INTORG/INTEND
// recognized and recorded but not yet consumed by any integer-aware solve
// path), RHS, RANGES, BOUNDS, ENDATA.
core::LpProblem ReadMps(const std::string& path);

}  // namespace inferno::io
