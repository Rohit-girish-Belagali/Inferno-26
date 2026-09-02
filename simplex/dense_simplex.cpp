#include "simplex/dense_simplex.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace inferno::simplex {

namespace {

using core::kInfinity;

enum class Status { kAtLower, kAtUpper, kFree, kBasic };

// Dense m x m Gauss-Jordan inversion with partial pivoting. Returns false if
// the matrix is numerically singular.
bool InvertDense(std::vector<std::vector<double>> a, std::vector<std::vector<double>>& inv,
                  double pivot_tol) {
  int m = static_cast<int>(a.size());
  inv.assign(m, std::vector<double>(m, 0.0));
  for (int i = 0; i < m; ++i) inv[i][i] = 1.0;

  for (int col = 0; col < m; ++col) {
    int pivot_row = col;
    double best = std::abs(a[col][col]);
    for (int r = col + 1; r < m; ++r) {
      if (std::abs(a[r][col]) > best) {
        best = std::abs(a[r][col]);
        pivot_row = r;
      }
    }
    if (best < pivot_tol) return false;
    if (pivot_row != col) {
      std::swap(a[pivot_row], a[col]);
      std::swap(inv[pivot_row], inv[col]);
    }
    double pivot = a[col][col];
    for (int c = 0; c < m; ++c) {
      a[col][c] /= pivot;
      inv[col][c] /= pivot;
    }
    for (int r = 0; r < m; ++r) {
      if (r == col) continue;
      double factor = a[r][col];
      if (factor == 0.0) continue;
      for (int c = 0; c < m; ++c) {
        a[r][c] -= factor * a[col][c];
        inv[r][c] -= factor * inv[col][c];
      }
    }
  }
  return true;
}

struct Workspace {
  int m = 0;  // rows = number of equality constraints (= problem.num_rows)
  int n = 0;  // total variables = num_cols + num_rows (structural + slack)
  std::vector<std::vector<double>> M;  // m x n dense constraint matrix, M w = 0
  std::vector<double> lo, hi, cost_phase2;
  std::vector<int> basis;       // size m, variable index basic in row k
  std::vector<Status> status;   // size n, meaningful for nonbasic entries
  std::vector<double> value;    // size n, current value of every variable
};

// alpha = Binv * column
std::vector<double> ApplyBinv(const std::vector<std::vector<double>>& binv,
                               const std::vector<double>& column) {
  int m = static_cast<int>(binv.size());
  std::vector<double> out(m, 0.0);
  for (int i = 0; i < m; ++i) {
    double s = 0.0;
    for (int k = 0; k < m; ++k) s += binv[i][k] * column[k];
    out[i] = s;
  }
  return out;
}

// One full pass of the primal simplex against a given cost vector, honoring
// per-variable status; recomputes x_B and Binv from scratch every
// iteration. Returns true if it terminated at optimality for `cost`, false
// on iteration-limit or singular basis (numerical_error_out set in that
// case).
bool RunSimplexPhase(Workspace& ws, const std::vector<double>& cost, bool is_phase1,
                      int max_iterations, const core::TolerancePolicy& tol,
                      bool& numerical_error, bool& hit_iteration_limit, bool& unbounded,
                      int& iterations_out) {
  numerical_error = false;
  hit_iteration_limit = false;
  unbounded = false;
  iterations_out = 0;

  for (int iter = 0; iter < max_iterations; ++iter) {
    iterations_out = iter + 1;
    std::vector<std::vector<double>> B(ws.m, std::vector<double>(ws.m));
    for (int k = 0; k < ws.m; ++k) {
      int var = ws.basis[k];
      for (int i = 0; i < ws.m; ++i) B[i][k] = ws.M[i][var];
    }
    std::vector<std::vector<double>> binv;
    if (!InvertDense(B, binv, tol.pivot)) {
      numerical_error = true;
      return false;
    }

    // x_B = -Binv * (M * x_with_basic_zeroed)
    std::vector<double> mx(ws.m, 0.0);
    for (int j = 0; j < ws.n; ++j) {
      if (ws.status[j] == Status::kBasic) continue;
      double xj = ws.value[j];
      if (xj == 0.0) continue;
      for (int i = 0; i < ws.m; ++i) mx[i] += ws.M[i][j] * xj;
    }
    std::vector<double> xb = ApplyBinv(binv, mx);
    for (int i = 0; i < ws.m; ++i) xb[i] = -xb[i];
    for (int k = 0; k < ws.m; ++k) ws.value[ws.basis[k]] = xb[k];

    // Effective phase-1 cost: -1 below lower, +1 above upper, else 0, for
    // basic variables only (nonbasic are feasible by construction).
    std::vector<double> cost_b(ws.m);
    for (int k = 0; k < ws.m; ++k) {
      int var = ws.basis[k];
      if (is_phase1) {
        double v = ws.value[var];
        if (std::isfinite(ws.lo[var]) && v < ws.lo[var] - tol.feasibility) {
          cost_b[k] = -1.0;
        } else if (std::isfinite(ws.hi[var]) && v > ws.hi[var] + tol.feasibility) {
          cost_b[k] = 1.0;
        } else {
          cost_b[k] = 0.0;
        }
      } else {
        cost_b[k] = cost[var];
      }
    }

    if (is_phase1) {
      double total_infeas = 0.0;
      for (int k = 0; k < ws.m; ++k) {
        int var = ws.basis[k];
        double v = ws.value[var];
        if (std::isfinite(ws.lo[var]) && v < ws.lo[var]) total_infeas += ws.lo[var] - v;
        if (std::isfinite(ws.hi[var]) && v > ws.hi[var]) total_infeas += v - ws.hi[var];
      }
      if (total_infeas <= tol.feasibility) return true;  // phase 1 complete, feasible
    }

    // y = Binv^T * cost_b : y[i] = sum_k binv[k][i] * cost_b[k]
    std::vector<double> y(ws.m, 0.0);
    for (int i = 0; i < ws.m; ++i) {
      double s = 0.0;
      for (int k = 0; k < ws.m; ++k) s += binv[k][i] * cost_b[k];
      y[i] = s;
    }

    // Pick entering variable via Bland's rule (smallest improving index).
    int entering = -1;
    int dir = 0;
    for (int j = 0; j < ws.n; ++j) {
      if (ws.status[j] == Status::kBasic) continue;
      double cj = is_phase1 ? 0.0 : cost[j];
      double reduced = cj;
      for (int i = 0; i < ws.m; ++i) reduced -= y[i] * ws.M[i][j];

      bool improving = false;
      int this_dir = 0;
      switch (ws.status[j]) {
        case Status::kAtLower:
          if (reduced < -tol.optimality) {
            improving = true;
            this_dir = +1;
          }
          break;
        case Status::kAtUpper:
          if (reduced > tol.optimality) {
            improving = true;
            this_dir = -1;
          }
          break;
        case Status::kFree:
          if (reduced < -tol.optimality) {
            improving = true;
            this_dir = +1;
          } else if (reduced > tol.optimality) {
            improving = true;
            this_dir = -1;
          }
          break;
        default:
          break;
      }
      if (improving) {
        entering = j;
        dir = this_dir;
        break;  // Bland: first improving index, no further scanning.
      }
    }

    if (entering == -1) return true;  // optimal for this phase/cost

    std::vector<double> alpha = ApplyBinv(binv, [&] {
      std::vector<double> col(ws.m);
      for (int i = 0; i < ws.m; ++i) col[i] = ws.M[i][entering];
      return col;
    }());

    double self_limit = kInfinity;
    if (std::isfinite(ws.lo[entering]) && std::isfinite(ws.hi[entering])) {
      self_limit = ws.hi[entering] - ws.lo[entering];
    }

    double best_t = self_limit;
    int leaving_row = -1;
    bool leaving_to_upper = false;

    for (int k = 0; k < ws.m; ++k) {
      double rate = -dir * alpha[k];
      if (std::abs(rate) < tol.pivot) continue;
      int var = ws.basis[k];
      double v = ws.value[var];
      double lo = ws.lo[var];
      double hi = ws.hi[var];

      bool infeasible_below = std::isfinite(lo) && v < lo - tol.feasibility;
      bool infeasible_above = std::isfinite(hi) && v > hi + tol.feasibility;

      double t_limit = kInfinity;
      bool candidate_to_upper = false;

      if (!infeasible_below && !infeasible_above) {
        if (rate > tol.pivot && std::isfinite(hi)) {
          t_limit = (hi - v) / rate;
          candidate_to_upper = true;
        } else if (rate < -tol.pivot && std::isfinite(lo)) {
          t_limit = (lo - v) / rate;
          candidate_to_upper = false;
        }
      } else if (infeasible_below) {
        if (rate > tol.pivot) {
          t_limit = (lo - v) / rate;
          candidate_to_upper = false;
        }
      } else {  // infeasible_above
        if (rate < -tol.pivot) {
          t_limit = (hi - v) / rate;
          candidate_to_upper = true;
        }
      }

      if (t_limit < -tol.feasibility) t_limit = 0.0;  // guard tiny negative noise
      if (t_limit < best_t) {
        best_t = t_limit;
        leaving_row = k;
        leaving_to_upper = candidate_to_upper;
      }
    }

    if (!std::isfinite(best_t)) {
      unbounded = true;
      return false;
    }

    double t = std::max(0.0, best_t);

    // Update all basic variable values along the step.
    for (int k = 0; k < ws.m; ++k) {
      ws.value[ws.basis[k]] += (-dir) * alpha[k] * t;
    }
    ws.value[entering] += dir * t;

    if (leaving_row == -1) {
      // Bound flip: entering variable itself stays nonbasic at the other bound.
      ws.status[entering] = (dir > 0) ? Status::kAtUpper : Status::kAtLower;
      ws.value[entering] = (dir > 0) ? ws.hi[entering] : ws.lo[entering];
    } else {
      int leaving_var = ws.basis[leaving_row];
      ws.status[leaving_var] = leaving_to_upper ? Status::kAtUpper : Status::kAtLower;
      ws.value[leaving_var] = leaving_to_upper ? ws.hi[leaving_var] : ws.lo[leaving_var];
      ws.status[entering] = Status::kBasic;
      ws.basis[leaving_row] = entering;
    }
  }

  hit_iteration_limit = true;
  return false;
}

}  // namespace

core::Solution SolveDense(const core::LpProblem& problem, int max_iterations,
                           const core::TolerancePolicy& tol) {
  core::Solution solution;

  Workspace ws;
  ws.m = problem.num_rows;
  ws.n = problem.num_cols + problem.num_rows;

  ws.lo.resize(ws.n);
  ws.hi.resize(ws.n);
  ws.cost_phase2.assign(ws.n, 0.0);
  ws.status.resize(ws.n);
  ws.value.assign(ws.n, 0.0);

  for (int j = 0; j < problem.num_cols; ++j) {
    ws.lo[j] = problem.col_lo[j];
    ws.hi[j] = problem.col_hi[j];
    ws.cost_phase2[j] = problem.obj[j];
  }
  for (int i = 0; i < problem.num_rows; ++i) {
    int slack = problem.num_cols + i;
    ws.lo[slack] = problem.row_lo[i];
    ws.hi[slack] = problem.row_hi[i];
    ws.cost_phase2[slack] = 0.0;
  }

  ws.M.assign(ws.m, std::vector<double>(ws.n, 0.0));
  for (int c = 0; c < problem.a.cols; ++c) {
    for (int p = problem.a.col_ptr[c]; p < problem.a.col_ptr[c + 1]; ++p) {
      ws.M[problem.a.row_idx[p]][c] = problem.a.values[p];
    }
  }
  for (int i = 0; i < problem.num_rows; ++i) {
    ws.M[i][problem.num_cols + i] = -1.0;
  }

  // Initial basis: all slacks. Nonbasic structural variables at their
  // nearest finite bound (lower preferred), or free (value 0) if unbounded
  // both ways.
  ws.basis.resize(ws.m);
  for (int i = 0; i < ws.m; ++i) ws.basis[i] = problem.num_cols + i;
  for (int i = 0; i < ws.m; ++i) ws.status[problem.num_cols + i] = Status::kBasic;

  for (int j = 0; j < problem.num_cols; ++j) {
    if (std::isfinite(ws.lo[j])) {
      ws.status[j] = Status::kAtLower;
      ws.value[j] = ws.lo[j];
    } else if (std::isfinite(ws.hi[j])) {
      ws.status[j] = Status::kAtUpper;
      ws.value[j] = ws.hi[j];
    } else {
      ws.status[j] = Status::kFree;
      ws.value[j] = 0.0;
    }
  }

  int iter_cap = max_iterations > 0 ? max_iterations : (200 * (ws.m + ws.n) + 2000);

  bool numerical_error = false, hit_limit = false, unbounded = false;
  int iters_used = 0;
  int phase1_iters = 0;
  std::vector<double> zero_cost(ws.n, 0.0);
  bool phase1_ok = RunSimplexPhase(ws, zero_cost, /*is_phase1=*/true, iter_cap, tol,
                                    numerical_error, hit_limit, unbounded, phase1_iters);
  iters_used += phase1_iters;

  if (numerical_error) {
    solution.status = core::SolveStatus::kNumericalError;
    return solution;
  }
  if (hit_limit) {
    solution.status = core::SolveStatus::kIterationLimit;
    return solution;
  }
  if (!phase1_ok) {
    // Should not happen (phase 1 objective is bounded below by 0), but be
    // defensive rather than silently reporting a wrong status.
    solution.status = core::SolveStatus::kNumericalError;
    return solution;
  }

  double total_infeas = 0.0;
  for (int j = 0; j < ws.n; ++j) {
    double v = ws.value[j];
    if (std::isfinite(ws.lo[j]) && v < ws.lo[j] - tol.feasibility) total_infeas += ws.lo[j] - v;
    if (std::isfinite(ws.hi[j]) && v > ws.hi[j] + tol.feasibility) total_infeas += v - ws.hi[j];
  }
  if (total_infeas > tol.checker_residual) {
    solution.status = core::SolveStatus::kInfeasible;
    return solution;
  }

  bool p2_numerical_error = false, p2_hit_limit = false, p2_unbounded = false;
  int phase2_iters = 0;
  bool phase2_ok = RunSimplexPhase(ws, ws.cost_phase2, /*is_phase1=*/false, iter_cap, tol,
                                    p2_numerical_error, p2_hit_limit, p2_unbounded, phase2_iters);
  iters_used += phase2_iters;

  if (p2_unbounded) {
    solution.status = core::SolveStatus::kUnbounded;
    return solution;
  }
  if (p2_numerical_error) {
    solution.status = core::SolveStatus::kNumericalError;
    return solution;
  }
  if (p2_hit_limit) {
    solution.status = core::SolveStatus::kIterationLimit;
    return solution;
  }
  (void)phase2_ok;

  // Final duals: recompute Binv one more time at the optimal basis.
  std::vector<std::vector<double>> B(ws.m, std::vector<double>(ws.m));
  for (int k = 0; k < ws.m; ++k) {
    int var = ws.basis[k];
    for (int i = 0; i < ws.m; ++i) B[i][k] = ws.M[i][var];
  }
  std::vector<std::vector<double>> binv;
  if (!InvertDense(B, binv, tol.pivot)) {
    solution.status = core::SolveStatus::kNumericalError;
    return solution;
  }
  std::vector<double> cost_b(ws.m);
  for (int k = 0; k < ws.m; ++k) cost_b[k] = ws.cost_phase2[ws.basis[k]];
  std::vector<double> y(ws.m, 0.0);
  for (int i = 0; i < ws.m; ++i) {
    double s = 0.0;
    for (int k = 0; k < ws.m; ++k) s += binv[k][i] * cost_b[k];
    y[i] = s;
  }

  solution.status = core::SolveStatus::kOptimal;
  solution.x.assign(problem.num_cols, 0.0);
  for (int j = 0; j < problem.num_cols; ++j) solution.x[j] = ws.value[j];

  solution.basis = ws.basis;

  solution.row_activity.assign(problem.num_rows, 0.0);
  for (int i = 0; i < problem.num_rows; ++i) solution.row_activity[i] = ws.value[problem.num_cols + i];

  solution.y = y;

  solution.reduced_cost.assign(problem.num_cols, 0.0);
  std::vector<double> aty(problem.num_cols, 0.0);
  problem.a.TransposeMultiplyAdd(y, aty);
  for (int j = 0; j < problem.num_cols; ++j) {
    solution.reduced_cost[j] = problem.obj[j] - aty[j];
  }

  double obj = problem.obj_offset;
  for (int j = 0; j < problem.num_cols; ++j) obj += problem.obj[j] * solution.x[j];
  solution.objective_value = obj;
  solution.iterations = iters_used;

  return solution;
}

}  // namespace inferno::simplex
