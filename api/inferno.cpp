#include "api/inferno.h"

#include <chrono>
#include <cmath>
#include <cstring>
#include <mutex>
#include <string>
#include <vector>

#include "checker/checker.hpp"
#include "checker/mip_checker.hpp"
#include "core/lp_problem.hpp"
#include "core/mip_problem.hpp"
#include "core/sparse.hpp"
#include "io/mps_reader.hpp"
#include "mip/branch_and_bound.hpp"
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

  // MILP state. is_integer stays empty for a pure LP.
  std::vector<char> is_integer;
  inferno::core::MipSolution mip;
  bool mip_solved = false;
  double solve_seconds = 0.0;

  // The event log, and the one piece of shared state in this ABI. A
  // caller observing a solve live reads it from another thread while the
  // search writes it, so both sides take this lock. It is deliberately
  // NOT a general concurrency story: only the two event accessors and the
  // solve callback touch it.
  mutable std::mutex event_mutex;
  std::vector<inferno_node_event> events;
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

// The inverse: the search legitimately produces IEEE infinities (no
// incumbent yet, no bound yet), and a C caller in another language may
// not be able to represent those. They cross the boundary as the same
// large finite sentinel every bound uses.
double ToApiBound(double v) {
  if (std::isnan(v)) return INFERNO_INFINITY;
  if (v >= INFERNO_INFINITY || std::isinf(v)) return v > 0 ? INFERNO_INFINITY : -INFERNO_INFINITY;
  if (v <= -INFERNO_INFINITY) return -INFERNO_INFINITY;
  return v;
}

}  // namespace

extern "C" {

inferno_status inferno_problem_set_integer(inferno_problem* p, const int* flags, int count) {
  if (p == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (flags == nullptr) {
    p->is_integer.clear();
    p->mip_solved = false;
    return INFERNO_STATUS_OK;
  }
  if (count != p->lp.num_cols) return INFERNO_STATUS_INVALID_ARGUMENT;
  try {
    p->is_integer.assign(static_cast<size_t>(count), 0);
    for (int j = 0; j < count; ++j) p->is_integer[static_cast<size_t>(j)] = flags[j] ? 1 : 0;
    p->mip_solved = false;
    return INFERNO_STATUS_OK;
  } catch (...) {
    return INFERNO_STATUS_INTERNAL_ERROR;
  }
}

inferno_status inferno_solve_mip(inferno_problem* p, double time_limit_seconds, int node_limit,
                                  double gap_tolerance, int presolve) {
  if (p == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!(time_limit_seconds > 0.0) || node_limit <= 0) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!(gap_tolerance >= 0.0)) return INFERNO_STATUS_INVALID_ARGUMENT;
  try {
    {
      std::lock_guard<std::mutex> lock(p->event_mutex);
      p->events.clear();
    }
    p->mip_solved = false;
    p->solved = false;

    inferno::core::MipProblem problem;
    problem.lp = p->lp;
    problem.is_integer = p->is_integer;
    if (problem.is_integer.empty()) {
      problem.is_integer.assign(static_cast<size_t>(p->lp.num_cols), 0);
    }

    inferno::mip::BranchAndBoundOptions opts;
    opts.time_limit_seconds = time_limit_seconds;
    opts.node_limit = node_limit;
    opts.gap_tolerance = gap_tolerance;
    opts.use_presolve = presolve != 0;
    opts.node_callback = [p](const inferno::mip::NodeEvent& ev) {
      inferno_node_event out;
      out.node_index = ev.node_index;
      out.depth = ev.depth;
      out.branch_var = ev.branch_var;
      out.branch_value = ev.branch_value;
      out.node_bound = ToApiBound(ev.node_bound);
      out.incumbent = ToApiBound(ev.incumbent);
      out.best_bound = ToApiBound(ev.best_bound);
      out.outcome = static_cast<int>(ev.outcome);
      out.elapsed_seconds = ev.elapsed_seconds;
      out.open_nodes = ev.open_nodes;
      std::lock_guard<std::mutex> lock(p->event_mutex);
      p->events.push_back(out);
    };

    const auto t0 = std::chrono::steady_clock::now();
    p->mip = inferno::mip::SolveMip(problem, opts);
    p->solve_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // Standing rule 2 again, at the MILP level: the independent MILP
    // checker decides whether this answer is acceptable, not the search
    // that produced it.
    p->checker_passed = false;
    if (p->mip.status == inferno::core::SolveStatus::kOptimal) {
      p->checker_passed = inferno::checker::VerifyMipSolution(problem, p->mip).passed;
    }

    // Mirror into the LP-shaped result so every existing accessor
    // (status, objective, solution vector) keeps working unchanged.
    p->solution = inferno::core::Solution{};
    p->solution.status = p->mip.status;
    p->solution.x = p->mip.x;
    p->solution.objective_value = p->mip.objective_value;
    p->solution.iterations = p->mip.nodes_explored;
    p->solved = true;
    p->mip_solved = true;
    return INFERNO_STATUS_OK;
  } catch (...) {
    return INFERNO_STATUS_INTERNAL_ERROR;
  }
}

inferno_status inferno_get_best_bound(const inferno_problem* p, double* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->mip_solved) return INFERNO_STATUS_NOT_SOLVED;
  *out = ToApiBound(p->mip.best_bound);
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_gap(const inferno_problem* p, double* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->mip_solved) return INFERNO_STATUS_NOT_SOLVED;
  *out = ToApiBound(p->mip.gap);
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_nodes_explored(const inferno_problem* p, int* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->mip_solved) return INFERNO_STATUS_NOT_SOLVED;
  *out = p->mip.nodes_explored;
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_proved_optimal(const inferno_problem* p, int* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->mip_solved) return INFERNO_STATUS_NOT_SOLVED;
  *out = p->mip.proved_optimal ? 1 : 0;
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_solve_seconds(const inferno_problem* p, double* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (!p->solved) return INFERNO_STATUS_NOT_SOLVED;
  *out = p->solve_seconds;
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_event_count(const inferno_problem* p, int* out) {
  if (p == nullptr || out == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  std::lock_guard<std::mutex> lock(p->event_mutex);
  *out = static_cast<int>(p->events.size());
  return INFERNO_STATUS_OK;
}

inferno_status inferno_get_events(const inferno_problem* p, inferno_node_event* out, int start,
                                   int count, int* written) {
  if (p == nullptr || out == nullptr || written == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  if (start < 0 || count < 0) return INFERNO_STATUS_INVALID_ARGUMENT;
  std::lock_guard<std::mutex> lock(p->event_mutex);
  const int have = static_cast<int>(p->events.size());
  int n = 0;
  while (n < count && start + n < have) {
    out[n] = p->events[static_cast<size_t>(start + n)];
    ++n;
  }
  *written = n;
  return INFERNO_STATUS_OK;
}

const char* inferno_node_outcome_string(int outcome) {
  switch (outcome) {
    case INFERNO_NODE_ROOT: return "root";
    case INFERNO_NODE_BRANCHED: return "branched";
    case INFERNO_NODE_INTEGER_FEASIBLE: return "integer";
    case INFERNO_NODE_INFEASIBLE: return "infeasible";
    case INFERNO_NODE_DOMINATED: return "dominated";
    case INFERNO_NODE_GAP_CUT: return "gap-cut";
    case INFERNO_NODE_RELAXATION_FAILED: return "unresolved";
  }
  return "unknown";
}

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
    p->mip_solved = false;
    p->is_integer.clear();
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
    p->mip_solved = false;
    p->is_integer.clear();
    return INFERNO_STATUS_OK;
  } catch (...) {
    return INFERNO_STATUS_INTERNAL_ERROR;
  }
}

inferno_status inferno_solve(inferno_problem* p, int presolve) {
  if (p == nullptr) return INFERNO_STATUS_INVALID_ARGUMENT;
  try {
    p->mip_solved = false;
    const auto t0 = std::chrono::steady_clock::now();
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
    p->solve_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
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

const char* inferno_version(void) { return "0.6.0"; }

}  // extern "C"
