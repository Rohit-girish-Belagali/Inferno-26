#include "io/mps_reader.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "core/lp_problem.hpp"
#include "core/sparse.hpp"

namespace inferno::io {

namespace {

using core::kInfinity;

enum class Section { kNone, kRows, kColumns, kRhs, kRanges, kBounds };

std::vector<std::string> Tokenize(const std::string& line) {
  std::vector<std::string> tokens;
  std::istringstream iss(line);
  std::string tok;
  while (iss >> tok) tokens.push_back(tok);
  return tokens;
}

bool IsCommentOrBlank(const std::string& line) {
  if (line.empty()) return true;
  std::size_t first = line.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return true;
  return line[first] == '*';
}

// True if the line starts in column 0 with a non-space character — that is
// the MPS convention for a section header (ROWS, COLUMNS, ...), as opposed
// to a data line which is indented.
bool IsSectionHeader(const std::string& line) {
  return !line.empty() && line[0] != ' ' && line[0] != '\t';
}

struct RowInfo {
  char type = 'N';
  int index = -1;  // index into LpProblem rows, -1 for the objective row
};

}  // namespace

core::LpProblem ReadMps(const std::string& path) {
  std::ifstream file(path);
  if (!file) {
    throw std::runtime_error("ReadMps: cannot open file: " + path);
  }

  core::LpProblem problem;
  Section section = Section::kNone;

  std::unordered_map<std::string, RowInfo> row_by_name;
  std::string objective_row_name;
  std::vector<std::string> row_order;  // constraint rows only, in file order

  std::unordered_map<std::string, int> col_index;
  std::vector<std::string> col_order;
  bool in_integer_block = false;
  std::vector<bool> col_is_integer;

  // Column entries collected as (col, row, value); converted to CSC once
  // the row count is final, i.e. after the ROWS section closes.
  struct Entry {
    int col;
    int row;  // -1 means objective
    double value;
  };
  std::vector<Entry> entries;

  std::unordered_map<std::string, double> rhs_by_row;   // row name -> rhs
  std::unordered_map<std::string, double> range_by_row;  // row name -> range value
  bool has_explicit_rhs_for_objective = false;
  double objective_rhs = 0.0;

  auto get_or_add_col = [&](const std::string& name) -> int {
    auto it = col_index.find(name);
    if (it != col_index.end()) return it->second;
    int idx = static_cast<int>(col_order.size());
    col_index.emplace(name, idx);
    col_order.push_back(name);
    col_is_integer.push_back(in_integer_block);
    return idx;
  };

  std::string line;
  int line_no = 0;
  while (std::getline(file, line)) {
    ++line_no;
    if (IsCommentOrBlank(line)) continue;

    if (IsSectionHeader(line)) {
      std::vector<std::string> header = Tokenize(line);
      const std::string& kw = header[0];
      if (kw == "NAME") {
        problem.name = header.size() > 1 ? header[1] : "";
        section = Section::kNone;
      } else if (kw == "ROWS") {
        section = Section::kRows;
      } else if (kw == "COLUMNS") {
        section = Section::kColumns;
      } else if (kw == "RHS") {
        section = Section::kRhs;
      } else if (kw == "RANGES") {
        section = Section::kRanges;
      } else if (kw == "BOUNDS") {
        section = Section::kBounds;
      } else if (kw == "ENDATA") {
        break;
      } else {
        // Unrecognized section (SOS, OBJSENSE, etc.) — skip its body.
        section = Section::kNone;
      }
      continue;
    }

    std::vector<std::string> tok = Tokenize(line);
    if (tok.empty()) continue;

    switch (section) {
      case Section::kRows: {
        if (tok.size() < 2) {
          throw std::runtime_error("ReadMps: malformed ROWS line " + std::to_string(line_no));
        }
        char type = tok[0][0];
        const std::string& name = tok[1];
        if (type == 'N') {
          if (objective_row_name.empty()) {
            objective_row_name = name;
            row_by_name[name] = RowInfo{type, -1};
          } else {
            // Extra free rows are accepted and ignored (not the objective).
            row_by_name[name] = RowInfo{type, -2};
          }
        } else {
          int idx = static_cast<int>(row_order.size());
          row_order.push_back(name);
          row_by_name[name] = RowInfo{type, idx};
        }
        break;
      }
      case Section::kColumns: {
        if (tok.size() >= 3 && tok[1] == "'MARKER'") {
          if (tok.size() >= 3 && tok[2] == "'INTORG'") {
            in_integer_block = true;
          } else if (tok.size() >= 3 && tok[2] == "'INTEND'") {
            in_integer_block = false;
          }
          break;
        }
        if (tok.size() < 3 || (tok.size() % 2) == 0) {
          throw std::runtime_error("ReadMps: malformed COLUMNS line " + std::to_string(line_no));
        }
        int col = get_or_add_col(tok[0]);
        for (std::size_t i = 1; i + 1 < tok.size(); i += 2) {
          const std::string& row_name = tok[i];
          double value = std::stod(tok[i + 1]);
          auto it = row_by_name.find(row_name);
          if (it == row_by_name.end()) {
            throw std::runtime_error("ReadMps: unknown row '" + row_name + "' at line " +
                                      std::to_string(line_no));
          }
          if (it->second.index == -1) {
            entries.push_back({col, -1, value});
          } else if (it->second.index >= 0) {
            entries.push_back({col, it->second.index, value});
          }
          // index == -2: extra free row, ignored.
        }
        break;
      }
      case Section::kRhs: {
        if (tok.size() < 3 || (tok.size() % 2) == 0) {
          throw std::runtime_error("ReadMps: malformed RHS line " + std::to_string(line_no));
        }
        for (std::size_t i = 1; i + 1 < tok.size(); i += 2) {
          const std::string& row_name = tok[i];
          double value = std::stod(tok[i + 1]);
          if (row_name == objective_row_name) {
            has_explicit_rhs_for_objective = true;
            objective_rhs = value;
          } else {
            rhs_by_row[row_name] = value;
          }
        }
        break;
      }
      case Section::kRanges: {
        if (tok.size() < 3 || (tok.size() % 2) == 0) {
          throw std::runtime_error("ReadMps: malformed RANGES line " + std::to_string(line_no));
        }
        for (std::size_t i = 1; i + 1 < tok.size(); i += 2) {
          range_by_row[tok[i]] = std::stod(tok[i + 1]);
        }
        break;
      }
      case Section::kBounds:
        // Bound values only matter once every column referenced in COLUMNS
        // has been seen, so BOUNDS is handled in a dedicated second pass
        // below rather than here.
        break;
      default:
        break;
    }
  }

  // Dedicated second pass for BOUNDS, now that col_index is fully populated.
  file.clear();
  file.seekg(0);
  section = Section::kNone;

  int num_cols = static_cast<int>(col_order.size());
  problem.col_lo.assign(num_cols, 0.0);
  problem.col_hi.assign(num_cols, kInfinity);

  line_no = 0;
  while (std::getline(file, line)) {
    ++line_no;
    if (IsCommentOrBlank(line)) continue;
    if (IsSectionHeader(line)) {
      std::vector<std::string> header = Tokenize(line);
      const std::string& kw = header[0];
      if (kw == "BOUNDS") {
        section = Section::kBounds;
      } else if (kw == "ENDATA") {
        break;
      } else {
        section = Section::kNone;
      }
      continue;
    }
    if (section != Section::kBounds) continue;

    std::vector<std::string> tok = Tokenize(line);
    if (tok.size() < 3) continue;
    const std::string& type = tok[0];
    const std::string& col_name = tok[2];
    auto it = col_index.find(col_name);
    if (it == col_index.end()) continue;
    int col = it->second;
    double value = tok.size() > 3 ? std::stod(tok[3]) : 0.0;

    if (type == "UP") {
      problem.col_hi[col] = value;
      if (value < 0.0 && problem.col_lo[col] == 0.0) problem.col_lo[col] = -kInfinity;
    } else if (type == "LO") {
      problem.col_lo[col] = value;
    } else if (type == "FX") {
      problem.col_lo[col] = value;
      problem.col_hi[col] = value;
    } else if (type == "FR") {
      problem.col_lo[col] = -kInfinity;
      problem.col_hi[col] = kInfinity;
    } else if (type == "MI") {
      problem.col_lo[col] = -kInfinity;
    } else if (type == "PL") {
      problem.col_hi[col] = kInfinity;
    } else if (type == "BV") {
      problem.col_lo[col] = 0.0;
      problem.col_hi[col] = 1.0;
    } else if (type == "UI") {
      problem.col_hi[col] = value;
    } else if (type == "LI") {
      problem.col_lo[col] = value;
    }
    // Unrecognized bound types are ignored rather than fatal, matching
    // common reader tolerance for vendor extensions.
  }

  // Assemble rows.
  problem.num_rows = static_cast<int>(row_order.size());
  problem.num_cols = num_cols;
  problem.row_lo.assign(problem.num_rows, -kInfinity);
  problem.row_hi.assign(problem.num_rows, kInfinity);
  problem.row_names = row_order;
  problem.col_names = col_order;

  for (int i = 0; i < problem.num_rows; ++i) {
    const std::string& name = row_order[i];
    char type = row_by_name.at(name).type;
    double rhs = 0.0;
    auto rhs_it = rhs_by_row.find(name);
    if (rhs_it != rhs_by_row.end()) rhs = rhs_it->second;

    switch (type) {
      case 'L':
        problem.row_lo[i] = -kInfinity;
        problem.row_hi[i] = rhs;
        break;
      case 'G':
        problem.row_lo[i] = rhs;
        problem.row_hi[i] = kInfinity;
        break;
      case 'E':
        problem.row_lo[i] = rhs;
        problem.row_hi[i] = rhs;
        break;
      default:
        throw std::runtime_error("ReadMps: unknown row type '" + std::string(1, type) +
                                  "' for row " + name);
    }

    auto range_it = range_by_row.find(name);
    if (range_it != range_by_row.end()) {
      double r = range_it->second;
      double abs_r = std::abs(r);
      switch (type) {
        case 'L':
          problem.row_lo[i] = rhs - abs_r;
          break;
        case 'G':
          problem.row_hi[i] = rhs + abs_r;
          break;
        case 'E':
          if (r >= 0.0) {
            problem.row_lo[i] = rhs;
            problem.row_hi[i] = rhs + abs_r;
          } else {
            problem.row_lo[i] = rhs - abs_r;
            problem.row_hi[i] = rhs;
          }
          break;
        default:
          break;
      }
    }
  }

  // Objective: minimize obj^T x + obj_offset. MPS convention is that a
  // nonzero RHS on the objective row is subtracted as a constant.
  problem.obj.assign(num_cols, 0.0);
  problem.obj_offset = has_explicit_rhs_for_objective ? -objective_rhs : 0.0;

  core::CscBuilder builder(problem.num_rows, num_cols);
  for (const auto& e : entries) {
    if (e.row == -1) {
      problem.obj[e.col] += e.value;
    } else {
      builder.AddEntry(e.col, e.row, e.value);
    }
  }
  problem.a = std::move(builder).Build();
  problem.a.Validate();

  if (!problem.IsFeasibleBounds()) {
    throw std::runtime_error("ReadMps: problem has crossed bounds (lo > hi) after parsing " +
                              path);
  }

  return problem;
}

}  // namespace inferno::io
