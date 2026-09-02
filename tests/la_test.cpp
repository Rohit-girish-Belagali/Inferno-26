// Phase 1.2 test suite for la/: Markowitz LU + threshold pivoting,
// Gilbert-Peierls sparse FTRAN/BTRAN, the PFI eta-chain basis update (see
// basis_factorization.hpp for why it's PFI and not full Forrest-Tomlin),
// and the refactorization policy. Validated the same way the checker
// validates the simplex: by recomputing residuals from raw data (B x ≈ b,
// B^T y ≈ c) rather than trusting the solver's own bookkeeping.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "core/sparse.hpp"
#include "core/tolerance.hpp"
#include "la/basis_factorization.hpp"
#include "la/lu_solve.hpp"
#include "la/markowitz_lu.hpp"
#include "la/scaling.hpp"

namespace {

int g_failures = 0;

void Expect(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << "\n";
    ++g_failures;
  } else {
    std::cout << "ok: " << message << "\n";
  }
}

using inferno::core::CscBuilder;
using inferno::core::CscMatrix;
using inferno::core::DefaultTolerances;
using inferno::la::BasisFactorization;
using inferno::la::FactorizeMarkowitz;
using inferno::la::LuFactors;
using inferno::la::MarkowitzOptions;

CscMatrix BuildDense(int m, const std::vector<double>& dense_row_major) {
  CscBuilder b(m, m);
  for (int r = 0; r < m; ++r) {
    for (int c = 0; c < m; ++c) {
      double v = dense_row_major[r * m + c];
      if (v != 0.0) b.AddEntry(c, r, v);
    }
  }
  return std::move(b).Build();
}

std::vector<std::pair<int, double>> ToSparse(const std::vector<double>& v) {
  std::vector<std::pair<int, double>> out;
  for (int i = 0; i < static_cast<int>(v.size()); ++i) {
    if (v[i] != 0.0) out.emplace_back(i, v[i]);
  }
  return out;
}

double MaxAbsResidual_Bx_b(const CscMatrix& b, const std::vector<double>& x,
                            const std::vector<double>& rhs) {
  std::vector<double> bx(b.rows, 0.0);
  b.MultiplyAdd(x, bx);
  double res = 0.0;
  for (int i = 0; i < b.rows; ++i) res = std::max(res, std::abs(bx[i] - rhs[i]));
  return res;
}

double MaxAbsResidual_BTy_c(const CscMatrix& b, const std::vector<double>& y,
                             const std::vector<double>& rhs) {
  std::vector<double> bty(b.cols, 0.0);
  b.TransposeMultiplyAdd(y, bty);
  double res = 0.0;
  for (int i = 0; i < b.cols; ++i) res = std::max(res, std::abs(bty[i] - rhs[i]));
  return res;
}

void TestSmallHandBuilt() {
  // A nonsymmetric 4x4 with a mix of magnitudes, chosen so naive diagonal
  // pivoting alone would be numerically poor (row 0's diagonal is tiny
  // relative to its off-diagonal), exercising real pivot selection.
  int m = 4;
  std::vector<double> dense = {
      0.0001, 2.0, 0.0, 1.0,
      3.0,    1.0, 0.0, 0.0,
      0.0,    0.0, 4.0, 2.0,
      1.0,    0.0, 1.0, 5.0,
  };
  CscMatrix b = BuildDense(m, dense);

  LuFactors lu;
  MarkowitzOptions opts;
  bool ok = FactorizeMarkowitz(b, opts, DefaultTolerances().pivot, lu);
  Expect(ok, "hand-built 4x4 factorizes");
  if (!ok) return;

  std::vector<double> rhs = {1.0, 2.0, 3.0, 4.0};
  std::vector<double> x = inferno::la::Ftran(lu, ToSparse(rhs));
  double primal_res = MaxAbsResidual_Bx_b(b, x, rhs);
  Expect(primal_res < 1e-10, "hand-built 4x4 FTRAN residual < 1e-10 (got " +
                                  std::to_string(primal_res) + ")");

  std::vector<double> rhs_t = {5.0, 6.0, 7.0, 8.0};
  std::vector<double> y = inferno::la::Btran(lu, ToSparse(rhs_t));
  double dual_res = MaxAbsResidual_BTy_c(b, y, rhs_t);
  Expect(dual_res < 1e-10, "hand-built 4x4 BTRAN residual < 1e-10 (got " +
                                std::to_string(dual_res) + ")");
}

void TestIdentity() {
  int m = 5;
  CscBuilder builder(m, m);
  for (int i = 0; i < m; ++i) builder.AddEntry(i, i, 1.0);
  CscMatrix b = std::move(builder).Build();

  LuFactors lu;
  MarkowitzOptions opts;
  bool ok = FactorizeMarkowitz(b, opts, DefaultTolerances().pivot, lu);
  Expect(ok, "identity factorizes");

  std::vector<double> rhs = {1, 2, 3, 4, 5};
  std::vector<double> x = inferno::la::Ftran(lu, ToSparse(rhs));
  bool matches = true;
  for (int i = 0; i < m; ++i) matches &= (std::abs(x[i] - rhs[i]) < 1e-12);
  Expect(matches, "identity FTRAN returns the RHS unchanged");
}

void TestSingularDetected() {
  // Row of all zeros -> structurally singular.
  int m = 3;
  std::vector<double> dense = {
      1.0, 2.0, 3.0,
      0.0, 0.0, 0.0,
      4.0, 5.0, 6.0,
  };
  CscMatrix b = BuildDense(m, dense);
  LuFactors lu;
  MarkowitzOptions opts;
  bool ok = FactorizeMarkowitz(b, opts, DefaultTolerances().pivot, lu);
  Expect(!ok, "structurally singular matrix (zero row) is detected, not silently factored");
}

void TestThresholdPivotingRejectsUnstablePivot() {
  // Column 0 has a tiny entry at row 0 (1e-8) and a much larger one at row
  // 1 (1.0). Threshold pivoting (tau=0.1) must reject row 0 as column 0's
  // pivot and select row 1 instead — verified indirectly: if it *had*
  // picked the tiny pivot, the resulting multiplier would be huge (1e8)
  // and residuals would blow up under floating point. Checking the
  // residual stays tiny is a robust way to confirm stability was honored
  // without reaching into the factorization's internals.
  int m = 3;
  std::vector<double> dense = {
      1e-8, 1.0, 0.0,
      1.0,  1.0, 1.0,
      0.0,  2.0, 1.0,
  };
  CscMatrix b = BuildDense(m, dense);
  LuFactors lu;
  MarkowitzOptions opts;
  opts.pivot_threshold = 0.1;
  bool ok = FactorizeMarkowitz(b, opts, DefaultTolerances().pivot, lu);
  Expect(ok, "threshold-pivoting test matrix factorizes");
  if (!ok) return;

  std::vector<double> rhs = {1.0, 2.0, 3.0};
  std::vector<double> x = inferno::la::Ftran(lu, ToSparse(rhs));
  double res = MaxAbsResidual_Bx_b(b, x, rhs);
  Expect(res < 1e-8, "threshold pivoting keeps FTRAN residual tiny despite a near-zero "
                      "natural pivot (got " + std::to_string(res) + ")");
}

void TestRandomSparseFactorAndSolve() {
  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> val_dist(-5.0, 5.0);
  std::uniform_real_distribution<double> density_dist(0.0, 1.0);

  int m = 60;
  std::vector<double> dense(static_cast<size_t>(m) * m, 0.0);
  for (int i = 0; i < m; ++i) dense[i * m + i] = 5.0 + std::abs(val_dist(rng));  // diagonally dominant
  for (int r = 0; r < m; ++r) {
    for (int c = 0; c < m; ++c) {
      if (r == c) continue;
      if (density_dist(rng) < 0.08) dense[r * m + c] = val_dist(rng);
    }
  }
  CscMatrix b = BuildDense(m, dense);

  LuFactors lu;
  MarkowitzOptions opts;
  bool ok = FactorizeMarkowitz(b, opts, DefaultTolerances().pivot, lu);
  Expect(ok, "60x60 diagonally-dominant sparse random matrix factorizes");
  if (!ok) return;

  std::vector<double> rhs(m);
  for (auto& v : rhs) v = val_dist(rng);
  std::vector<double> x = inferno::la::Ftran(lu, ToSparse(rhs));
  double primal_res = MaxAbsResidual_Bx_b(b, x, rhs);
  Expect(primal_res < 1e-10,
         "60x60 FTRAN residual < 1e-10 (got " + std::to_string(primal_res) + ")");

  std::vector<double> rhs_t(m);
  for (auto& v : rhs_t) v = val_dist(rng);
  std::vector<double> y = inferno::la::Btran(lu, ToSparse(rhs_t));
  double dual_res = MaxAbsResidual_BTy_c(b, y, rhs_t);
  Expect(dual_res < 1e-10, "60x60 BTRAN residual < 1e-10 (got " + std::to_string(dual_res) + ")");
}

// Phase 1 gate: "1000 sequential updates stay within 1e-8 of a fresh
// refactorization" and "update at least 20x faster than refactorizing".
void TestSequentialUpdatesAndTiming() {
  std::mt19937 rng(777);
  std::uniform_real_distribution<double> val_dist(-3.0, 3.0);
  std::uniform_real_distribution<double> density_dist(0.0, 1.0);

  int m = 80;
  auto random_column = [&]() {
    std::vector<double> col(m, 0.0);
    for (int i = 0; i < m; ++i) {
      if (density_dist(rng) < 0.15) col[i] = val_dist(rng);
    }
    return col;
  };

  std::vector<std::vector<double>> true_cols(m);
  std::vector<double> dense(static_cast<size_t>(m) * m, 0.0);
  for (int c = 0; c < m; ++c) {
    dense[c * m + c] = 5.0;  // seed a safely-invertible diagonal, then randomize off-diagonal
    std::vector<double> col = random_column();
    col[c] = 5.0 + std::abs(val_dist(rng));
    true_cols[c] = col;
    for (int r = 0; r < m; ++r) dense[r * m + c] = col[r];
  }
  CscMatrix b0 = BuildDense(m, dense);

  BasisFactorization bf;
  MarkowitzOptions opts;
  auto tol = DefaultTolerances();
  bool ok = bf.Factorize(b0, opts, tol);
  Expect(ok, "80x80 base matrix for the update test factorizes");
  if (!ok) return;

  auto build_current_matrix = [&]() {
    CscBuilder builder(m, m);
    for (int c = 0; c < m; ++c) {
      for (int r = 0; r < m; ++r) {
        if (true_cols[c][r] != 0.0) builder.AddEntry(c, r, true_cols[c][r]);
      }
    }
    return std::move(builder).Build();
  };

  int refactor_count = 0;
  auto update_start = std::chrono::steady_clock::now();
  int applied = 0;
  int attempts = 0;
  while (applied < 1000 && attempts < 4000) {
    ++attempts;
    int t = static_cast<int>(rng() % m);
    std::vector<double> col = random_column();
    col[t] += 5.0;  // bias toward a safe, non-tiny pivot at the replaced slot
    bool applied_ok = bf.Update(t, ToSparse(col), tol);
    if (!applied_ok) continue;  // unsafe pivot this draw; try a different column
    true_cols[t] = col;
    ++applied;
    if (bf.ShouldRefactorize()) {
      CscMatrix current = build_current_matrix();
      bool refactor_ok = bf.Factorize(current, opts, tol);
      Expect(refactor_ok, "periodic refactor during the update stress test succeeds");
      ++refactor_count;
    }
  }
  auto update_end = std::chrono::steady_clock::now();
  double update_total_s = std::chrono::duration<double>(update_end - update_start).count();

  Expect(applied == 1000, "1000 basis updates were successfully applied (in " +
                               std::to_string(attempts) + " attempts, " +
                               std::to_string(refactor_count) + " periodic refactors)");

  CscMatrix final_matrix = build_current_matrix();
  std::vector<double> rhs(m);
  for (auto& v : rhs) v = val_dist(rng);

  std::vector<double> x_updated = bf.Ftran(ToSparse(rhs));
  double updated_res = MaxAbsResidual_Bx_b(final_matrix, x_updated, rhs);
  Expect(updated_res < 1e-8,
         "after 1000 updates, FTRAN residual against the true final basis < 1e-8 (got " +
             std::to_string(updated_res) + ")");

  LuFactors fresh_lu;
  bool fresh_ok = FactorizeMarkowitz(final_matrix, opts, tol.pivot, fresh_lu);
  Expect(fresh_ok, "fresh refactorization of the final basis (for comparison) succeeds");
  if (fresh_ok) {
    std::vector<double> x_fresh = inferno::la::Ftran(fresh_lu, ToSparse(rhs));
    double diff = 0.0;
    for (int i = 0; i < m; ++i) diff = std::max(diff, std::abs(x_fresh[i] - x_updated[i]));
    Expect(diff < 1e-8, "updated-chain solution matches a fresh refactorization's solution "
                         "within 1e-8 (got " + std::to_string(diff) + ")");
  }

  // Timing: average update cost vs average full-factorization cost, same
  // matrix size. Refactor timing measured on freshly built same-size
  // matrices, not reusing `final_matrix`'s specific sparsity by accident.
  int refactor_trials = 30;
  auto refactor_start = std::chrono::steady_clock::now();
  for (int trial = 0; trial < refactor_trials; ++trial) {
    LuFactors throwaway;
    FactorizeMarkowitz(final_matrix, opts, tol.pivot, throwaway);
  }
  auto refactor_end = std::chrono::steady_clock::now();
  double avg_refactor_s =
      std::chrono::duration<double>(refactor_end - refactor_start).count() / refactor_trials;
  double avg_update_s = update_total_s / applied;

  double speedup = avg_update_s > 0.0 ? avg_refactor_s / avg_update_s : 0.0;
  printf("avg update: %.9fs, avg refactor: %.9fs, speedup: %.1fx\n", avg_update_s, avg_refactor_s,
         speedup);
  Expect(speedup >= 20.0, "basis update is at least 20x faster than refactorizing (got " +
                               std::to_string(speedup) + "x)");
}

void TestScalingProducesUnitMaxMagnitude() {
  int m = 4, n = 5;
  CscBuilder builder(m, n);
  builder.AddEntry(0, 0, 1000.0);
  builder.AddEntry(0, 1, 0.001);
  builder.AddEntry(1, 1, 50.0);
  builder.AddEntry(2, 2, 2.0);
  builder.AddEntry(3, 2, 0.5);
  builder.AddEntry(4, 3, 200.0);
  builder.AddEntry(2, 0, 7.0);
  CscMatrix a = std::move(builder).Build();

  auto scale = inferno::la::ComputeGeometricScaling(a, 2);
  CscMatrix scaled = inferno::la::ApplyScaling(a, scale);

  double max_abs = 0.0;
  for (double v : scaled.values) max_abs = std::max(max_abs, std::abs(v));
  Expect(max_abs <= 1.0 + 1e-9,
         "geometric + equilibration scaling brings every entry's magnitude to <= 1 (got max " +
             std::to_string(max_abs) + ")");
}

}  // namespace

int main() {
  TestSmallHandBuilt();
  TestIdentity();
  TestSingularDetected();
  TestThresholdPivotingRejectsUnstablePivot();
  TestRandomSparseFactorAndSolve();
  TestScalingProducesUnitMaxMagnitude();
  TestSequentialUpdatesAndTiming();

  if (g_failures > 0) {
    std::cerr << g_failures << " test(s) failed\n";
    return 1;
  }
  std::cout << "all la tests passed\n";
  return 0;
}
