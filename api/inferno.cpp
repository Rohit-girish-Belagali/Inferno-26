#include "api/inferno.h"

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "checker/checker.hpp"
#include "core/lp_problem.hpp"
#include "core/sparse.hpp"
#include "io/mps_reader.hpp"
#include "presolve/presolve.hpp"
#include "simplex/revised_simplex.hpp"

// The opaque handle. Keeping the C++ types entirely inside this
// translation unit is the point of the whole file: no C++ type ever
// appears in the header, so the ABI does not move when the standard
// library or compiler does.
struct inferno_problem {
  inferno::core::LpProblem lp;
  inferno::core::Solution solution;
  bool solved = false;
  bool checker_passed = false;
};

namespace {

// INFERNO_INFINITY is a large finite double rather than HUGE_VAL, because
// a caller in another language may not have a clean way to express IEEE
// infinity. Anything at or beyond it means unbounded.
double FromApiBound(double v) {
  if (v >= INFERNO_INFINITY) return inferno::core::kInfinity;
  if (v <= -INFERNO_INFINITY) return -inferno::core::kInfinity;
  return v;
}

}  // namespace

extern "C" {

inferno_problem* inferno_problem_create(void) {
  // Nothing here may throw into the caller, and operator new can.
  try {
    return new inferno_problem();
  } catch (...) {
    return nullptr;
  }
}

void inferno_problem_destroy(inferno_problem* p) { delete p; }

inferno_status inferno_problem_load_mps(inferno_problem* p, const char* path) {
  if (p == nullptr || path == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  try {
    inferno::core::LpProblem parsed = inferno::io::ReadMps(path);
    p->lp = std::move(parsed);  // only on success, so a failed read leaves p untouched
    p->solved = false;
    return INFERNO_STATUS_OK;
  } catch (const std::exception&) {
    return INFERNO_STATUS_READ_FAILED;
  } catch (...) {
    return INFERNO_STATUS_INTERNAL_ERROR;
  }
}

inferno_status inferno_problem_set(inferno_problem* p, int num_rows, int num_cols,
                                    const int* col_ptr, const int* row_idx,
                                    const double* values, const double* obj,
                                    const double* col_lo, const double* col_hi,
                                    const double* row_lo, const double* row_hi) {
  if (p == nullptr || col_ptr == nullptr || obj == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (col_lo == nullptr || col_hi == nullptr || row_lo == nullptr || row_hi == nullptr) {
    return INFERNO_STATUS_INVALID_ARGUMENT;
  }
  if (num_rows < 0 || num_cols < 0) return INFERNO_STATUS_INVALID_ARGUMENT;
  const int nnz = num_cols > 0 ? col_ptr[num_cols] : 0;
  if (nnz < 0) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (nnz > 0 && (row_idx == nullptr || values == nullptr)) return INFERNO_STATUS_INVALID_ARGUMENT;

  try {
    // Validate before mutating anything: a caller that passes a malformed
    // matrix should get an error and keep whatever problem it had, not a
    // half-written one.
    for (int c = 0; c < num_cols; ++c) {
      if (col_ptr[c] > col_ptr[c + 1] || col_ptr[c] < 0) return INFERNO_STATUS_INVALID_ARGUMENT;
    }
    for (int k = 0; k < nnz; ++k) {
      if (row_idx[k] < 0 || row_idx[k] >= num_rows) return INFERNO_STATUS_INVALID_ARGUMENT;
    }

    inferno::core::LpProblem lp;
    lp.name = "API";
    lp.num_rows = num_rows;
    lp.num_cols = num_cols;
    inferno::core::CscBuilder b(num_rows, num_cols);
    for (int c = 0; c < num_cols; ++c) {
      for (int k = col_ptr[c]; k < col_ptr[c + 1]; ++k) b.AddEntry(c, row_idx[k], values[k]);
    }
    lp.a = std::move(b).Build();
    lp.obj.assign(obj, obj + num_cols);
    lp.col_lo.resize(num_cols);
    lp.col_hi.resize(num_cols);
    for (int j = 0; j < num_cols; ++j) {
      lp.col_lo[j] = FromApiBound(col_lo[j]);
      lp.col_hi[j] = FromApiBound(col_hi[j]);
    }
    lp.row_lo.resize(num_rows);
    lp.row_hi.resize(num_rows);
    for (int i = 0; i < num_rows; ++i) {
      lp.row_lo[i] = FromApiBound(row_lo[i]);
      lp.row_hi[i] = FromApiBound(row_hi[i]);
    }
    lp.col_names.assign(num_cols, std::string());
    lp.row_names.assign(num_rows, std::string());

    p->lp = std::move(lp);
    p->solved = false;
    return INFERNO_STATUS_OK;
  } catch (...) {
    return INFERNO_STATUS_INTERNAL_ERROR;
  }
}

inferno_status inferno_solve(inferno_problem* p, int presolve) {
  if (p == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  try {
    if (presolve) {
      auto pre = inferno::presolve::Presolve(p->lp);
      if (pre.infeasible) {
        p->solution = inferno::core::Solution{};
        p->solution.status = inferno::core::SolveStatus::kInfeasible;
      } else {
        inferno::core::Solution reduced = inferno::simplex::SolveRevised(pre.reduced);
        p->solution = inferno::presolve::Postsolve(p->lp, pre, reduced);
      }
    } else {
      p->solution = inferno::simplex::SolveRevised(p->lp);
    }

    // Standing rule 2: the checker runs on every solve, in every mode.
    // Exposing it through the ABI keeps that true for embedders too --
    // otherwise the discipline would stop at the language boundary.
    p->checker_passed = false;
    if (p->solution.status == inferno::core::SolveStatus::kOptimal) {
      p->checker_passed = inferno::checker::VerifySolution(p->lp, p->solution).passed;
    }
    p->solved = true;
    return INFERNO_STATUS_OK;
  } catch (...) {
    return INFERNO_STATUS_INTERNAL_ERROR;
  }
}

inferno_status inferno_get_solve_status(const inferno_problem* p, inferno_solve_status* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->solved) return INFERNO_STATUS_NOT_SOLVED;
  switch (p->solution.status) {
    case inferno::core::SolveStatus::kOptimal: *out = INFERNO_SOLVE_OPTIMAL; break;
    case inferno::core::SolveStatus::kInfeasible: *out = INFERNO_SOLVE_INFEASIBLE; break;
    case inferno::core::SolveStatus::kUnbounded: *out = INFERNO_SOLVE_UNBOUNDED; break;
    case inferno::core::SolveStatus::kIterationLimit: *out = INFERNO_SOLVE_ITERATION_LIMIT; break;
    default: *out = INFERNO_SOLVE_NUMERICAL_ERROR; break;
  }
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_objective(const inferno_problem* p, double* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->solved) return INFERNO_STATUS_NOT_SOLVED;
  *out = p->solution.objective_value;
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_iterations(const inferno_problem* p, int* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->solved) return INFERNO_STATUS_NOT_SOLVED;
  *out = p->solution.iterations;
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_checker_passed(const inferno_problem* p, int* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->solved) return INFERNO_STATUS_NOT_SOLVED;
  *out = p->checker_passed ? 1 : 0;
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_solution(const inferno_problem* p, double* out, int count) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->solved) return INFERNO_STATUS_NOT_SOLVED;
  if (count != p->lp.num_cols) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (static_cast<int>(p->solution.x.size()) != count) return INFERNO_STATUS_NOT_SOLVED;
  std::memcpy(out, p->solution.x.data(), sizeof(double) * static_cast<size_t>(count));
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_duals(const inferno_problem* p, double* out, int count) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->solved) return INFERNO_STATUS_NOT_SOLVED;
  if (count != p->lp.num_rows) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (static_cast<int>(p->solution.y.size()) != count) return INFERNO_STATUS_NOT_SOLVED;
  std::memcpy(out, p->solution.y.data(), sizeof(double) * static_cast<size_t>(count));
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_num_rows(const inferno_problem* p, int* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  *out = p->lp.num_rows;
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_num_cols(const inferno_problem* p, int* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  *out = p->lp.num_cols;
  return INFERNO_STATUS_OK;
}

const char* inferno_status_string(inferno_status s) {
  switch (s) {
    case INFERNO_STATUS_OK: return "ok";
    case INFERNO_STATUS_INVALID_ARGUMENT: return "invalid argument";
    case INFERNO_STATUS_READ_FAILED: return "could not read or parse the problem file";
    case INFERNO_STATUS_NOT_SOLVED: return "solve has not been run on this problem yet";
    case INFERNO_STATUS_INTERNAL_ERROR: return "internal error";
  }
  return "unknown status";
}

const char* inferno_version(void) { return "0.5.0"; }

}  // extern "C"
