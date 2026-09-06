#include "qp/admm.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

#include "core/sparse.hpp"
#include "la/lu_solve.hpp"
#include "la/markowitz_lu.hpp"

namespace inferno::qp {

namespace {

using core::kInfinity;

double Clamp(double v, double lo, double hi) {
  if (std::isfinite(lo) && v < lo) return lo;
  if (std::isfinite(hi) && v > hi) return hi;
  return v;
}

// y += P x, where only P's upper triangle is stored. The strict upper part
// must be applied twice — once as itself and once transposed — because the
// stored half stands in for both.
void SymMultiply(const core::CscMatrix& p_upper, const std::vector<double>& x,
                  std::vector<double>& y) {
  for (int c = 0; c < p_upper.cols; ++c) {
    for (int idx = p_upper.col_ptr[c]; idx < p_upper.col_ptr[c + 1]; ++idx) {
      int r = p_upper.row_idx[idx];
      double v = p_upper.values[idx];
      y[r] += v * x[c];
      if (r != c) y[c] += v * x[r];
    }
  }
}

// Assembles the KKT matrix [P + sigma I, A'; A, -1/rho I] as a general
// (n+m) square sparse matrix. Both triangles are written out: the LU is a
// general factorization, not a symmetric one, so it needs the whole thing.
core::CscMatrix BuildKkt(const core::QpProblem& prob, double sigma, double rho) {
  const int n = prob.num_cols, m = prob.num_rows;

  // Entries are ACCUMULATED per position before being handed to
  // CscBuilder, which matters more than it looks: CscBuilder appends
  // whatever it is given and does not merge duplicates, and the LU then
  // builds its working matrix by map ASSIGNMENT — so a second entry at an
  // occupied position silently overwrites the first rather than adding to
  // it. P's diagonal and the sigma regularization land on exactly the same
  // position, so adding them as two entries replaced a diagonal of 1 with
  // one of 1e-6, leaving a near-singular block that drove the solution
  // straight to its bound. Merging here is the fix.
  std::vector<std::map<int, double>> cols(n + m);
  auto acc = [&](int col, int row, double v) { cols[col][row] += v; };

  for (int c = 0; c < n; ++c) {
    for (int idx = prob.p_upper.col_ptr[c]; idx < prob.p_upper.col_ptr[c + 1]; ++idx) {
      int r = prob.p_upper.row_idx[idx];
      double v = prob.p_upper.values[idx];
      acc(c, r, v);
      if (r != c) acc(r, c, v);  // mirror, since only the upper triangle is stored
    }
  }
  for (int c = 0; c < n; ++c) acc(c, c, sigma);

  // A in the lower-left block, A' in the upper-right.
  for (int c = 0; c < n; ++c) {
    for (int idx = prob.a.col_ptr[c]; idx < prob.a.col_ptr[c + 1]; ++idx) {
      int r = prob.a.row_idx[idx];
      double v = prob.a.values[idx];
      acc(c, n + r, v);
      acc(n + r, c, v);
    }
  }
  for (int i = 0; i < m; ++i) acc(n + i, n + i, -1.0 / rho);

  core::CscBuilder b(n + m, n + m);
  for (int c = 0; c < n + m; ++c) {
    for (const auto& [r, v] : cols[c]) {
      if (v != 0.0) b.AddEntry(c, r, v);
    }
  }
  return std::move(b).Build();
}

}  // namespace

core::QpSolution SolveQpAdmm(const core::QpProblem& problem, const AdmmOptions& opts,
                              const core::TolerancePolicy& tol) {
  core::QpSolution out;
  const int n = problem.num_cols, m = problem.num_rows;
  out.x.assign(n, 0.0);
  out.y.assign(m, 0.0);
  if (n == 0) { out.status = core::SolveStatus::kOptimal; return out; }

  double rho = opts.rho;
  la::LuFactors lu;
  auto factorize = [&](double r) {
    core::CscMatrix kkt = BuildKkt(problem, opts.sigma, r);
    return la::FactorizeMarkowitz(kkt, la::MarkowitzOptions{}, tol.pivot, lu);
  };
  if (!factorize(rho)) { out.status = core::SolveStatus::kNumericalError; return out; }

  std::vector<double> x(n, 0.0), z(m, 0.0), y(m, 0.0);
  std::vector<double> ax(m, 0.0), px(n, 0.0), aty(n, 0.0);

  int iter = 0;
  for (; iter < opts.max_iterations; ++iter) {
    // --- x-step: solve the KKT system. ---
    std::vector<std::pair<int, double>> rhs;
    rhs.reserve(n + m);
    for (int j = 0; j < n; ++j) {
      double v = opts.sigma * x[j] - problem.q[j];
      if (v != 0.0) rhs.emplace_back(j, v);
    }
    for (int i = 0; i < m; ++i) {
      double v = z[i] - y[i] / rho;
      if (v != 0.0) rhs.emplace_back(n + i, v);
    }
    std::vector<double> sol = la::Ftran(lu, rhs);

    std::vector<double> x_t(n), z_t(m);
    for (int j = 0; j < n; ++j) x_t[j] = sol[j];
    // nu comes back in the lower block; z~ = z + (nu - y)/rho.
    for (int i = 0; i < m; ++i) z_t[i] = z[i] + (sol[n + i] - y[i]) / rho;

    // --- Relaxation, projection, dual update. ---
    std::vector<double> z_prev = z;
    for (int j = 0; j < n; ++j) x[j] = opts.alpha * x_t[j] + (1.0 - opts.alpha) * x[j];
    for (int i = 0; i < m; ++i) {
      double relaxed = opts.alpha * z_t[i] + (1.0 - opts.alpha) * z_prev[i];
      // The y^k/rho term belongs inside the PROJECTION but not inside the
      // dual update. Carrying it into both makes the two terms cancel
      // exactly, so y never decays toward zero on an inactive constraint —
      // it holds whatever value it picked up and drags x away with it.
      // (Measured: with the extra term this diverged to x = -2e5.)
      z[i] = Clamp(relaxed + y[i] / rho, problem.row_lo[i], problem.row_hi[i]);
      y[i] += rho * (relaxed - z[i]);
    }

    if ((iter + 1) % opts.check_every != 0) continue;

    // --- KKT residuals. ---
    std::fill(ax.begin(), ax.end(), 0.0);
    problem.a.MultiplyAdd(x, ax);
    double primal = 0.0, ax_scale = 1.0;
    for (int i = 0; i < m; ++i) {
      double viol = 0.0;
      if (std::isfinite(problem.row_lo[i])) viol = std::max(viol, problem.row_lo[i] - ax[i]);
      if (std::isfinite(problem.row_hi[i])) viol = std::max(viol, ax[i] - problem.row_hi[i]);
      primal = std::max(primal, viol);
      ax_scale = std::max(ax_scale, std::abs(ax[i]));
    }

    std::fill(px.begin(), px.end(), 0.0);
    SymMultiply(problem.p_upper, x, px);
    std::fill(aty.begin(), aty.end(), 0.0);
    problem.a.TransposeMultiplyAdd(y, aty);
    double dual = 0.0, d_scale = 1.0;
    for (int j = 0; j < n; ++j) {
      double g = px[j] + problem.q[j] + aty[j];
      dual = std::max(dual, std::abs(g));
      d_scale = std::max({d_scale, std::abs(px[j]), std::abs(problem.q[j]), std::abs(aty[j])});
    }

    double rel_p = primal / ax_scale, rel_d = dual / d_scale;
    if (rel_p <= tol.checker_residual && rel_d <= tol.checker_residual) break;

    // --- Adaptive rho. When one residual runs far ahead of the other the
    // penalty is mismatched to the problem's scaling; rebalance and
    // refactorize. Only on a large discrepancy, since refactorizing is the
    // one genuinely expensive thing in this loop. ---
    if (opts.adaptive_rho && rel_p > 0.0 && rel_d > 0.0) {
      double ratio = std::sqrt(rel_p / rel_d);
      if (ratio > 5.0 || ratio < 0.2) {
        double new_rho = std::min(std::max(rho * ratio, 1e-6), 1e6);
        if (factorize(new_rho)) rho = new_rho;
      }
    }
  }

  out.iterations = iter;
  out.x = x;
  out.y = y;

  std::fill(ax.begin(), ax.end(), 0.0);
  problem.a.MultiplyAdd(x, ax);
  for (int i = 0; i < m; ++i) {
    double viol = 0.0;
    if (std::isfinite(problem.row_lo[i])) viol = std::max(viol, problem.row_lo[i] - ax[i]);
    if (std::isfinite(problem.row_hi[i])) viol = std::max(viol, ax[i] - problem.row_hi[i]);
    out.primal_residual = std::max(out.primal_residual, viol);
  }
  std::fill(px.begin(), px.end(), 0.0);
  SymMultiply(problem.p_upper, x, px);
  std::fill(aty.begin(), aty.end(), 0.0);
  problem.a.TransposeMultiplyAdd(y, aty);
  for (int j = 0; j < n; ++j) {
    out.dual_residual = std::max(out.dual_residual, std::abs(px[j] + problem.q[j] + aty[j]));
  }
  // Complementarity: a nonzero dual demands its constraint be active on the
  // matching side.
  for (int i = 0; i < m; ++i) {
    if (y[i] > tol.optimality && std::isfinite(problem.row_hi[i])) {
      out.complementarity = std::max(out.complementarity, std::abs(ax[i] - problem.row_hi[i]));
    } else if (y[i] < -tol.optimality && std::isfinite(problem.row_lo[i])) {
      out.complementarity = std::max(out.complementarity, std::abs(ax[i] - problem.row_lo[i]));
    }
  }

  double obj = problem.obj_offset;
  for (int j = 0; j < n; ++j) obj += problem.q[j] * x[j];
  obj += 0.5 * [&] { double t = 0.0; for (int j = 0; j < n; ++j) t += x[j] * px[j]; return t; }();
  out.objective_value = obj;

  bool ok = out.primal_residual <= tol.checker_residual &&
             out.dual_residual <= tol.checker_residual &&
             out.complementarity <= tol.checker_residual;
  out.status = ok ? core::SolveStatus::kOptimal : core::SolveStatus::kIterationLimit;
  return out;
}

}  // namespace inferno::qp
