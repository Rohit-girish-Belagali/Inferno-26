#include "simplex/dual_simplex.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/sparse.hpp"
#include "la/basis_factorization.hpp"
#include "la/markowitz_lu.hpp"
#include "la/scaling.hpp"

namespace inferno::simplex {

namespace {

using core::kInfinity;

enum class Status { kAtLower, kAtUpper, kFree, kBasic };

std::vector<std::pair<int, double>> ColumnOf(const core::LpProblem& problem, int var) {
  std::vector<std::pair<int, double>> col;
  if (var < problem.num_cols) {
    for (int p = problem.a.col_ptr[var]; p < problem.a.col_ptr[var + 1]; ++p) {
      col.emplace_back(problem.a.row_idx[p], problem.a.values[p]);
    }
  } else {
    col.emplace_back(var - problem.num_cols, -1.0);
  }
  return col;
}

core::CscMatrix BuildBasisMatrix(const core::LpProblem& problem, const std::vector<int>& basis) {
  int m = problem.num_rows;
  core::CscBuilder builder(m, m);
  for (int slot = 0; slot < m; ++slot) {
    for (const auto& [row, val] : ColumnOf(problem, basis[slot])) builder.AddEntry(slot, row, val);
  }
  return std::move(builder).Build();
}

struct Workspace {
  int m = 0, n = 0;
  std::vector<double> lo, hi, cost;
  std::vector<int> basis, basis_slot_of;
  std::vector<Status> status;
  std::vector<double> value;
  la::BasisFactorization bf;
};

void RecomputeBasicValues(Workspace& ws, const core::LpProblem& problem) {
  std::vector<double> mx(ws.m, 0.0);
  for (int j = 0; j < ws.n; ++j) {
    if (ws.basis_slot_of[j] != -1) continue;
    double xj = ws.value[j];
    if (xj == 0.0) continue;
    for (const auto& [row, val] : ColumnOf(problem, j)) mx[row] += val * xj;
  }
  std::vector<std::pair<int, double>> rhs_sparse;
  for (int row = 0; row < ws.m; ++row) {
    if (mx[row] != 0.0) rhs_sparse.emplace_back(row, -mx[row]);
  }
  std::vector<double> xb = ws.bf.Ftran(rhs_sparse);
  for (int slot = 0; slot < ws.m; ++slot) ws.value[ws.basis[slot]] = xb[slot];
}

double MaxBoundViolation(const Workspace& ws) {
  double worst = 0.0;
  for (int j = 0; j < ws.n; ++j) {
    double v = ws.value[j];
    if (std::isfinite(ws.lo[j])) worst = std::max(worst, ws.lo[j] - v);
    if (std::isfinite(ws.hi[j])) worst = std::max(worst, v - ws.hi[j]);
  }
  return worst;
}

// Solves `problem` (already scaled or not — the caller decides) with no
// further transformation. SolveDual() below applies scaling around this,
// same split as revised_simplex.cpp.
core::Solution SolveDualCore(const core::LpProblem& problem, int max_iterations,
                              const core::TolerancePolicy& tol) {
  core::Solution solution;

  Workspace ws;
  ws.m = problem.num_rows;
  ws.n = problem.num_cols + problem.num_rows;
  ws.lo.resize(ws.n);
  ws.hi.resize(ws.n);
  ws.cost.assign(ws.n, 0.0);
  ws.status.resize(ws.n);
  ws.value.assign(ws.n, 0.0);
  ws.basis.resize(ws.m);
  ws.basis_slot_of.assign(ws.n, -1);

  for (int j = 0; j < problem.num_cols; ++j) {
    ws.lo[j] = problem.col_lo[j];
    ws.hi[j] = problem.col_hi[j];
    ws.cost[j] = problem.obj[j];
  }
  for (int i = 0; i < problem.num_rows; ++i) {
    int slack = problem.num_cols + i;
    ws.lo[slack] = problem.row_lo[i];
    ws.hi[slack] = problem.row_hi[i];
  }

  for (int i = 0; i < ws.m; ++i) {
    int slack = problem.num_cols + i;
    ws.basis[i] = slack;
    ws.basis_slot_of[slack] = i;
    ws.status[slack] = Status::kBasic;
  }

  // Trivial dual-feasible start: with the all-slack basis, y = 0, so every
  // nonbasic structural variable's reduced cost is just its own objective
  // coefficient. Place it at whichever bound that sign requires; fail
  // (documented scope limit, see dual_simplex.hpp) if that bound doesn't
  // exist.
  for (int j = 0; j < problem.num_cols; ++j) {
    bool lo_finite = std::isfinite(ws.lo[j]);
    bool hi_finite = std::isfinite(ws.hi[j]);
    if (!lo_finite && !hi_finite) {
      if (std::abs(ws.cost[j]) <= tol.optimality) {
        ws.status[j] = Status::kFree;
        ws.value[j] = 0.0;
      } else {
        solution.status = core::SolveStatus::kNumericalError;
        return solution;
      }
    } else if (ws.cost[j] >= -tol.optimality) {
      if (!lo_finite) {
        solution.status = core::SolveStatus::kNumericalError;
        return solution;
      }
      ws.status[j] = Status::kAtLower;
      ws.value[j] = ws.lo[j];
    } else {
      if (!hi_finite) {
        solution.status = core::SolveStatus::kNumericalError;
        return solution;
      }
      ws.status[j] = Status::kAtUpper;
      ws.value[j] = ws.hi[j];
    }
  }

  std::vector<double> x_nonbasic(problem.num_cols);
  for (int j = 0; j < problem.num_cols; ++j) x_nonbasic[j] = ws.value[j];
  std::vector<double> ax(ws.m, 0.0);
  problem.a.MultiplyAdd(x_nonbasic, ax);
  for (int i = 0; i < ws.m; ++i) ws.value[problem.num_cols + i] = ax[i];

  la::MarkowitzOptions opts;
  core::CscMatrix b0 = BuildBasisMatrix(problem, ws.basis);
  if (!ws.bf.Factorize(b0, opts, tol)) {
    solution.status = core::SolveStatus::kNumericalError;
    return solution;
  }

  int iter_cap = max_iterations > 0 ? max_iterations : (200 * (ws.m + ws.n) + 2000);
  int iters_used = 0;
  int degenerate_streak = 0;
  constexpr int kBlandThreshold = 50;

  for (int iter = 0; iter < iter_cap; ++iter) {
    iters_used = iter + 1;
    // Anti-cycling fallback, mirroring revised_simplex.cpp's primal
    // pricing: after a run of degenerate (zero-length) pivots, switch both
    // selections to smallest-index-first (Bland's rule) instead of
    // largest-infeasibility / min-ratio, which guarantees termination.
    bool use_bland = degenerate_streak >= kBlandThreshold;

    int leaving_slot = -1;
    double worst = tol.feasibility;
    bool leaving_below = false;
    for (int slot = 0; slot < ws.m; ++slot) {
      int var = ws.basis[slot];
      double v = ws.value[var];
      bool below = std::isfinite(ws.lo[var]) && ws.lo[var] - v > tol.feasibility;
      bool above = std::isfinite(ws.hi[var]) && v - ws.hi[var] > tol.feasibility;
      if (!below && !above) continue;
      double viol = below ? (ws.lo[var] - v) : (v - ws.hi[var]);
      if (use_bland) {
        leaving_slot = slot;
        leaving_below = below;
        break;
      }
      if (viol > worst) {
        worst = viol;
        leaving_slot = slot;
        leaving_below = below;
      }
    }

    if (leaving_slot == -1) {
      // Primal-feasible while dual feasibility was maintained throughout
      // -> optimal. Ground first so this conclusion reflects the true
      // state, not incremental drift (same defensive pattern as
      // revised_simplex.cpp's phase-2 optimality check).
      RecomputeBasicValues(ws, problem);
      if (MaxBoundViolation(ws) > tol.checker_residual) {
        solution.status = core::SolveStatus::kNumericalError;
        return solution;
      }
      break;
    }

    std::vector<double> cost_b(ws.m);
    for (int slot = 0; slot < ws.m; ++slot) cost_b[slot] = ws.cost[ws.basis[slot]];
    std::vector<std::pair<int, double>> cost_b_sparse;
    for (int slot = 0; slot < ws.m; ++slot) {
      if (cost_b[slot] != 0.0) cost_b_sparse.emplace_back(slot, cost_b[slot]);
    }
    std::vector<double> y = ws.bf.Btran(cost_b_sparse);
    std::vector<double> rho = ws.bf.Btran({{leaving_slot, 1.0}});

    int entering = -1;
    int dir = 0;
    double best_ratio = 0.0;
    double best_alpha_abs = -1.0;

    for (int j = 0; j < ws.n; ++j) {
      if (ws.basis_slot_of[j] != -1) continue;
      if (ws.status[j] == Status::kFree) continue;  // not handled — see header scope note

      double alpha_rj = 0.0;
      for (const auto& [row, val] : ColumnOf(problem, j)) alpha_rj += rho[row] * val;
      if (std::abs(alpha_rj) < tol.pivot) continue;

      bool eligible;
      int this_dir;
      if (ws.status[j] == Status::kAtLower) {
        this_dir = +1;
        eligible = leaving_below ? (alpha_rj < -tol.pivot) : (alpha_rj > tol.pivot);
      } else {
        this_dir = -1;
        eligible = leaving_below ? (alpha_rj > tol.pivot) : (alpha_rj < -tol.pivot);
      }
      if (!eligible) continue;

      if (use_bland) {
        entering = j;
        dir = this_dir;
        break;
      }

      double reduced = ws.cost[j];
      for (const auto& [row, val] : ColumnOf(problem, j)) reduced -= y[row] * val;
      double ratio = std::abs(reduced / alpha_rj);

      if (entering == -1 || ratio < best_ratio - tol.optimality ||
          (std::abs(ratio - best_ratio) <= tol.optimality && std::abs(alpha_rj) > best_alpha_abs)) {
        best_ratio = ratio;
        best_alpha_abs = std::abs(alpha_rj);
        entering = j;
        dir = this_dir;
      }
    }

    if (entering == -1) {
      // Dual unbounded in the direction needed to fix this row's
      // infeasibility -> by LP duality, the primal is infeasible.
      solution.status = core::SolveStatus::kInfeasible;
      return solution;
    }

    std::vector<std::pair<int, double>> entering_col = ColumnOf(problem, entering);
    std::vector<double> alpha = ws.bf.Ftran(entering_col);

    int leaving_var = ws.basis[leaving_slot];
    double target = leaving_below ? ws.lo[leaving_var] : ws.hi[leaving_var];
    double rate = -dir * alpha[leaving_slot];
    if (std::abs(rate) < tol.pivot) {
      solution.status = core::SolveStatus::kNumericalError;
      return solution;
    }
    double t = (target - ws.value[leaving_var]) / rate;
    if (t < -tol.feasibility) {
      solution.status = core::SolveStatus::kNumericalError;
      return solution;
    }
    t = std::max(0.0, t);
    if (t < tol.feasibility) ++degenerate_streak; else degenerate_streak = 0;

    for (int slot = 0; slot < ws.m; ++slot) ws.value[ws.basis[slot]] += (-dir) * alpha[slot] * t;
    ws.value[entering] += dir * t;

    ws.status[leaving_var] = leaving_below ? Status::kAtLower : Status::kAtUpper;
    ws.value[leaving_var] = target;
    ws.basis_slot_of[leaving_var] = -1;
    ws.status[entering] = Status::kBasic;
    ws.basis_slot_of[entering] = leaving_slot;
    ws.basis[leaving_slot] = entering;

    bool update_ok = ws.bf.Update(leaving_slot, entering_col, tol);
    if (!update_ok || ws.bf.ShouldRefactorize()) {
      core::CscMatrix current = BuildBasisMatrix(problem, ws.basis);
      if (!ws.bf.Factorize(current, la::MarkowitzOptions{}, tol)) {
        solution.status = core::SolveStatus::kNumericalError;
        return solution;
      }
      RecomputeBasicValues(ws, problem);
    }

    if (iter == iter_cap - 1) {
      solution.status = core::SolveStatus::kIterationLimit;
      return solution;
    }
  }

  std::vector<double> cost_b(ws.m);
  for (int slot = 0; slot < ws.m; ++slot) cost_b[slot] = ws.cost[ws.basis[slot]];
  std::vector<std::pair<int, double>> cost_b_sparse;
  for (int slot = 0; slot < ws.m; ++slot) {
    if (cost_b[slot] != 0.0) cost_b_sparse.emplace_back(slot, cost_b[slot]);
  }
  std::vector<double> y = ws.bf.Btran(cost_b_sparse);

  solution.status = core::SolveStatus::kOptimal;
  solution.x.assign(problem.num_cols, 0.0);
  for (int j = 0; j < problem.num_cols; ++j) solution.x[j] = ws.value[j];

  solution.row_activity.assign(problem.num_rows, 0.0);
  for (int i = 0; i < problem.num_rows; ++i) {
    solution.row_activity[i] = ws.value[problem.num_cols + i];
  }

  solution.y = y;
  solution.basis = ws.basis;

  solution.reduced_cost.assign(problem.num_cols, 0.0);
  std::vector<double> aty(problem.num_cols, 0.0);
  problem.a.TransposeMultiplyAdd(y, aty);
  for (int j = 0; j < problem.num_cols; ++j) solution.reduced_cost[j] = problem.obj[j] - aty[j];

  double obj = problem.obj_offset;
  for (int j = 0; j < problem.num_cols; ++j) obj += problem.obj[j] * solution.x[j];
  solution.objective_value = obj;
  solution.iterations = iters_used;

  return solution;
}

}  // namespace

core::Solution SolveDual(const core::LpProblem& problem, int max_iterations,
                          const core::TolerancePolicy& tol) {
  la::ScaleFactors scale = la::ComputeGeometricScaling(problem.a);

  core::LpProblem scaled = problem;
  scaled.a = la::ApplyScaling(problem.a, scale);
  for (int i = 0; i < problem.num_rows; ++i) {
    scaled.row_lo[i] =
        std::isfinite(problem.row_lo[i]) ? problem.row_lo[i] * scale.row_scale[i] : problem.row_lo[i];
    scaled.row_hi[i] =
        std::isfinite(problem.row_hi[i]) ? problem.row_hi[i] * scale.row_scale[i] : problem.row_hi[i];
  }
  for (int j = 0; j < problem.num_cols; ++j) {
    scaled.col_lo[j] =
        std::isfinite(problem.col_lo[j]) ? problem.col_lo[j] / scale.col_scale[j] : problem.col_lo[j];
    scaled.col_hi[j] =
        std::isfinite(problem.col_hi[j]) ? problem.col_hi[j] / scale.col_scale[j] : problem.col_hi[j];
    scaled.obj[j] = problem.obj[j] * scale.col_scale[j];
  }

  core::Solution scaled_solution = SolveDualCore(scaled, max_iterations, tol);

  core::Solution solution;
  solution.status = scaled_solution.status;
  solution.iterations = scaled_solution.iterations;
  solution.basis = scaled_solution.basis;

  if (scaled_solution.status != core::SolveStatus::kOptimal) return solution;

  solution.x.assign(problem.num_cols, 0.0);
  for (int j = 0; j < problem.num_cols; ++j) solution.x[j] = scaled_solution.x[j] * scale.col_scale[j];

  solution.row_activity.assign(problem.num_rows, 0.0);
  for (int i = 0; i < problem.num_rows; ++i) {
    solution.row_activity[i] = scaled_solution.row_activity[i] / scale.row_scale[i];
  }

  solution.y.assign(problem.num_rows, 0.0);
  for (int i = 0; i < problem.num_rows; ++i) solution.y[i] = scaled_solution.y[i] * scale.row_scale[i];

  solution.reduced_cost.assign(problem.num_cols, 0.0);
  for (int j = 0; j < problem.num_cols; ++j) {
    solution.reduced_cost[j] = scaled_solution.reduced_cost[j] / scale.col_scale[j];
  }

  double obj = problem.obj_offset;
  for (int j = 0; j < problem.num_cols; ++j) obj += problem.obj[j] * solution.x[j];
  solution.objective_value = obj;

  return solution;
}

}  // namespace inferno::simplex
