#pragma once

#include <string>

#include "core/lp_problem.hpp"

namespace inferno::io {

// Reads an MPS file into a core::LpProblem. Tries free-form (whitespace-
// delimited fields) first, since most files including most of the Netlib
// set parse that way; on failure, retries in strict fixed-column mode
// (standard field layout: columns 2-3, 5-12, 15-22, 25-36, 40-47, 50-61),
// which is needed for the minority of older Netlib instances whose names
// contain embedded spaces or that leave a leading field legitimately blank
// on a continuation line — both of which whitespace tokenization corrupts.
// Throws std::runtime_error, with both modes' errors included, if neither
// mode parses the file.
//
// Supported sections: NAME, ROWS, COLUMNS (with MARKER INTORG/INTEND
// recognized and recorded but not yet consumed by any integer-aware solve
// path), RHS, RANGES, BOUNDS, ENDATA.
core::LpProblem ReadMps(const std::string& path);

}  // namespace inferno::io
