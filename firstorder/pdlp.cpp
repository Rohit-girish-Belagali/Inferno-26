#include "firstorder/pdlp.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/sparse.hpp"
#include "la/scaling.hpp"

namespace inferno::firstorder {

namespace {

using core::kInfinity;

double Clamp(double v, double lo, double hi) {
  if (std::isfinite(lo) && v < lo) return lo;
  if (std::isfinite(hi) && v > hi) return hi;
  return v;
}

// ||A||_2, estimated by power iteration on A'A. PDHG's step sizes are only
// valid below 1/||A||^2, so an UNDER-estimate would break convergence
// outright while an over-estimate merely slows it — hence the deliberate
// safety margin applied by the caller.
double EstimateNorm(const core::CscMatrix& a, int rows, int cols) {
  std::vector<double> v(cols, 1.0 / std::sqrt(static_cast<double>(std::max(cols, 1))));
  std::vector<double> av(rows), atav(cols);
  double norm = 1.0;
  for (int it = 0; it < 100; ++it) {
    std::fill(av.begin(), av.end(), 0.0);
    a.MultiplyAdd(v, av);
    std::fill(atav.begin(), atav.end(), 0.0);
    a.TransposeMultiplyAdd(av, atav);
    double n = 0.0;
    for (double x : atav) n += x * x;
    n = std::sqrt(n);
    if (n <= 0.0) return 1.0;
    double prev = norm;
    norm = std::sqrt(n);
    for (int j = 0; j < cols; ++j) v[j] = atav[j] / n;
    if (std::abs(norm - prev) <= 1e-10 * std::max(1.0, norm)) break;
  }
  return std::max(norm, 1e-12);
}

}  // namespace

PdlpResult SolvePdlpCore(const core::LpProblem& problem, const PdlpOptions& opts,
                          const core::TolerancePolicy& tol) {
  PdlpResult result;
  const int m = problem.num_rows, n = problem.num_cols;
  const core::CscMatrix& a = problem.a;

  std::vector<double> x(n, 0.0), y(m, 0.0);
  for (int j = 0; j < n; ++j) x[j] = Clamp(0.0, problem.col_lo[j], problem.col_hi[j]);

  // Running average of the iterates. PDHG's ergodic (averaged) sequence has
  // the better convergence guarantee, and restarting from it is what makes
  // the restart scheme work at all.
  std::vector<double> x_sum(n, 0.0), y_sum(m, 0.0);
  double sum_weight = 0.0;

  const double a_norm = EstimateNorm(a, m, n);
  // tau * sigma * ||A||^2 < 1 is required; 0.9 keeps a margin against the
  // power iteration having under-estimated the norm.
  double primal_weight = 1.0;  // ratio sigma/tau, rebalanced below
  double step = 0.9 / a_norm;

  std::vector<double> ax(m), aty(n), x_new(n), y_new(m), x_bar(n), scratch_m(m);
  x_bar = x;

  double best_kkt = kInfinity;
  std::vector<double> best_x = x, best_y = y;

  int iter = 0;
  for (; iter < opts.max_iterations; ++iter) {
    const double tau = step / primal_weight;   // primal step
    const double sigma = step * primal_weight; // dual step

    // --- Dual ascent: y+ = prox_{sigma h*}(y + sigma A xbar), and by
    // Moreau prox_{sigma h*}(v) = v - sigma proj_C(v / sigma). ---
    std::fill(ax.begin(), ax.end(), 0.0);
    a.MultiplyAdd(x_bar, ax);
    for (int i = 0; i < m; ++i) {
      double v = y[i] + sigma * ax[i];
      double proj = Clamp(v / sigma, problem.row_lo[i], problem.row_hi[i]);
      y_new[i] = v - sigma * proj;
    }

    // --- Primal descent: x+ = proj_X(x - tau (c + A' y+)). ---
    std::fill(aty.begin(), aty.end(), 0.0);
    a.TransposeMultiplyAdd(y_new, aty);
    for (int j = 0; j < n; ++j) {
      double g = problem.obj[j] + aty[j];
      x_new[j] = Clamp(x[j] - tau * g, problem.col_lo[j], problem.col_hi[j]);
    }

    for (int j = 0; j < n; ++j) x_bar[j] = 2.0 * x_new[j] - x[j];
    x.swap(x_new);
    y.swap(y_new);

    for (int j = 0; j < n; ++j) x_sum[j] += x[j];
    for (int i = 0; i < m; ++i) y_sum[i] += y[i];
    sum_weight += 1.0;

    if ((iter + 1) % opts.check_every != 0) continue;

    // --- Convergence test on the AVERAGED iterate. ---
    std::vector<double> xa(n), ya(m);
    for (int j = 0; j < n; ++j) xa[j] = x_sum[j] / sum_weight;
    for (int i = 0; i < m; ++i) ya[i] = y_sum[i] / sum_weight;

    std::fill(scratch_m.begin(), scratch_m.end(), 0.0);
    a.MultiplyAdd(xa, scratch_m);
    double primal_res = 0.0, row_scale = 1.0;
    for (int i = 0; i < m; ++i) {
      double viol = 0.0;
      if (std::isfinite(problem.row_lo[i])) viol = std::max(viol, problem.row_lo[i] - scratch_m[i]);
      if (std::isfinite(problem.row_hi[i])) viol = std::max(viol, scratch_m[i] - problem.row_hi[i]);
      primal_res = std::max(primal_res, viol);
      row_scale = std::max(row_scale, std::abs(scratch_m[i]));
    }

    // Dual residual: the reduced cost c + A'y must be sign-consistent with
    // where x sits in its box. Any part that is not is dual infeasibility.
    std::vector<double> atya(n, 0.0);
    a.TransposeMultiplyAdd(ya, atya);
    double dual_res = 0.0, obj_scale = 1.0;
    double dual_obj = 0.0;
    for (int j = 0; j < n; ++j) {
      double z = problem.obj[j] + atya[j];
      obj_scale = std::max(obj_scale, std::abs(problem.obj[j]));
      // z must be >= 0 where x can still rise, <= 0 where it can still fall.
      bool can_rise = !std::isfinite(problem.col_hi[j]) || xa[j] < problem.col_hi[j] - tol.feasibility;
      bool can_fall = !std::isfinite(problem.col_lo[j]) || xa[j] > problem.col_lo[j] + tol.feasibility;
      if (can_rise && z < 0.0) dual_res = std::max(dual_res, -z);
      if (can_fall && z > 0.0) dual_res = std::max(dual_res, z);
      // Dual objective contribution from the column box.
      if (z > 0.0 && std::isfinite(problem.col_lo[j])) dual_obj += z * problem.col_lo[j];
      else if (z < 0.0 && std::isfinite(problem.col_hi[j])) dual_obj += z * problem.col_hi[j];
    }
    // -h*(y), where h* is the support function of the row box: its
    // supremum sits at row_hi where y is positive and row_lo where y is
    // negative. If the bound it wants is infinite the support function is
    // unbounded, meaning this y is dual infeasible and the gap it implies
    // is meaningless — flag it rather than quietly reporting a small gap.
    // The sign test needs a tolerance, not an exact comparison. The prox
    // step keeps each y_i on its feasible side by construction, but these
    // are AVERAGED iterates, so a component whose true value is zero shows
    // up as numerical dust of arbitrary sign — and treating a -1e-18 on a
    // row with no lower bound as "wants an infinite bound" declared the
    // whole dual infeasible and reported an infinite gap on instances that
    // were otherwise converging perfectly well.
    bool dual_feasible = true;
    const double y_zero = 1e-11;
    for (int i = 0; i < m; ++i) {
      if (ya[i] > y_zero) {
        if (!std::isfinite(problem.row_hi[i])) { dual_feasible = false; break; }
        dual_obj -= ya[i] * problem.row_hi[i];
      } else if (ya[i] < -y_zero) {
        if (!std::isfinite(problem.row_lo[i])) { dual_feasible = false; break; }
        dual_obj -= ya[i] * problem.row_lo[i];
      }
    }

    double primal_obj = problem.obj_offset;
    for (int j = 0; j < n; ++j) primal_obj += problem.obj[j] * xa[j];
    double gap = dual_feasible ? std::abs(primal_obj - (dual_obj + problem.obj_offset)) : kInfinity;

    double rel_p = primal_res / row_scale;
    double rel_d = dual_res / obj_scale;
    double rel_g = gap / std::max(1.0, std::abs(primal_obj));
    double kkt = std::max({rel_p, rel_d, rel_g});

    if (kkt < best_kkt) {
      best_kkt = kkt;
      best_x = xa;
      best_y = ya;
    }

    result.relative_primal_residual = rel_p;
    result.relative_dual_residual = rel_d;
    result.relative_duality_gap = rel_g;

    // Terminate against the TIGHTER of the two tolerances, not the
    // checker's own bound. The KKT residuals measured here are relative
    // and computed on the averaged iterate, while the checker's are
    // absolute and recomputed independently from the raw problem — so
    // stopping the moment the relative measure crosses 1e-6 produced
    // solutions the checker then rejected (afiro stopped at 1151
    // iterations and failed; left to reach this tighter target it
    // converges to 1e-15 and passes). Claiming kOptimal is only worth
    // anything if the independent check agrees.
    if (kkt <= tol.feasibility) break;

    // --- Adaptive restart: collapse to the average and start fresh. The
    // averaged iterate has the stronger guarantee, so restarting the
    // non-averaged sequence from it is what converts PDHG's slow ergodic
    // rate into something practical. ---
    if (opts.use_restarts) {
      x = xa;
      y = ya;
      x_bar = xa;
      std::fill(x_sum.begin(), x_sum.end(), 0.0);
      std::fill(y_sum.begin(), y_sum.end(), 0.0);
      sum_weight = 0.0;
    }

    // --- Primal weight balancing. If the primal residual dominates the
    // dual one the primal steps are too timid relative to the dual, and
    // vice versa; nudge the ratio toward whichever side is lagging. Moved
    // geometrically and clamped, because reacting sharply to one noisy
    // measurement destabilizes the iteration. ---
    if (opts.use_primal_weight_balancing && rel_p > 0.0 && rel_d > 0.0) {
      double ratio = std::sqrt(rel_p / rel_d);
      primal_weight *= std::pow(ratio, 0.2);
      primal_weight = std::min(std::max(primal_weight, 1e-4), 1e4);
    }
  }

  result.iterations = iter;

  core::Solution& sol = result.solution;
  sol.x = best_x;
  // Sign convention. The saddle-point form used above puts the dual term
  // as +y'(Ax), so its reduced cost is c + A'y. Everything else in this
  // project — the checker included — uses z = c - A'y. Those differ by
  // the sign of y, so flip it here rather than leaving two conventions
  // live in the same codebase for a reader to trip over.
  sol.y.assign(m, 0.0);
  for (int i = 0; i < m; ++i) sol.y[i] = -best_y[i];
  sol.row_activity.assign(m, 0.0);
  a.MultiplyAdd(sol.x, sol.row_activity);
  sol.reduced_cost.assign(n, 0.0);
  std::vector<double> aty_final(n, 0.0);
  a.TransposeMultiplyAdd(sol.y, aty_final);
  // This project's sign convention (see checker/): z = c - A'y.
  for (int j = 0; j < n; ++j) sol.reduced_cost[j] = problem.obj[j] - aty_final[j];
  double obj = problem.obj_offset;
  for (int j = 0; j < n; ++j) obj += problem.obj[j] * sol.x[j];
  sol.objective_value = obj;
  sol.iterations = iter;
  sol.status = (best_kkt <= tol.feasibility) ? core::SolveStatus::kOptimal
                                              : core::SolveStatus::kIterationLimit;
  return result;
}



// Preconditioning wrapper. A first-order method's convergence rate depends
// directly on the conditioning of A — far more sharply than simplex's
// does, since there is no factorization absorbing the scale differences.
// The plan lists Ruiz preconditioning for exactly this reason; this reuses
// the geometric-mean + equilibration scaling already built in `la/` for
// the simplex path rather than adding a second, near-identical scaling
// implementation beside it.
//
// Measured on the instances that failed without it: `blend` went from not
// converging at all (relative primal residual 2.4e-2 after 200k
// iterations) to solving, and `share2b` from actively diverging (2.0e+01)
// to converging. The core iteration was already correct before this —
// `afiro` and `sc50a` reached the simplex objective to 1e-7 unscaled — so
// this is about conditioning, not correctness.
PdlpResult SolvePdlp(const core::LpProblem& problem, const PdlpOptions& opts,
                      const core::TolerancePolicy& tol) {
  la::ScaleFactors scale = la::ComputeGeometricScaling(problem.a);

  core::LpProblem scaled = problem;
  scaled.a = la::ApplyScaling(problem.a, scale);
  for (int i = 0; i < problem.num_rows; ++i) {
    scaled.row_lo[i] = std::isfinite(problem.row_lo[i]) ? problem.row_lo[i] * scale.row_scale[i]
                                                         : problem.row_lo[i];
    scaled.row_hi[i] = std::isfinite(problem.row_hi[i]) ? problem.row_hi[i] * scale.row_scale[i]
                                                         : problem.row_hi[i];
  }
  for (int j = 0; j < problem.num_cols; ++j) {
    scaled.col_lo[j] = std::isfinite(problem.col_lo[j]) ? problem.col_lo[j] / scale.col_scale[j]
                                                         : problem.col_lo[j];
    scaled.col_hi[j] = std::isfinite(problem.col_hi[j]) ? problem.col_hi[j] / scale.col_scale[j]
                                                         : problem.col_hi[j];
    scaled.obj[j] = problem.obj[j] * scale.col_scale[j];
  }

  PdlpResult r = SolvePdlpCore(scaled, opts, tol);

  core::Solution& s = r.solution;
  for (int j = 0; j < problem.num_cols; ++j) s.x[j] *= scale.col_scale[j];
  for (int i = 0; i < problem.num_rows; ++i) s.y[i] *= scale.row_scale[i];
  s.row_activity.assign(problem.num_rows, 0.0);
  problem.a.MultiplyAdd(s.x, s.row_activity);
  s.reduced_cost.assign(problem.num_cols, 0.0);
  std::vector<double> aty(problem.num_cols, 0.0);
  problem.a.TransposeMultiplyAdd(s.y, aty);
  for (int j = 0; j < problem.num_cols; ++j) s.reduced_cost[j] = problem.obj[j] - aty[j];
  double obj = problem.obj_offset;
  for (int j = 0; j < problem.num_cols; ++j) obj += problem.obj[j] * s.x[j];
  s.objective_value = obj;
  return r;
}

}  // namespace inferno::firstorder
