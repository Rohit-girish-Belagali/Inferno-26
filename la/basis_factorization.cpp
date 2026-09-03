#include "la/basis_factorization.hpp"

#include <algorithm>
#include <cmath>

#include "la/lu_solve.hpp"

namespace inferno::la {

bool BasisFactorization::Factorize(const core::CscMatrix& b, const MarkowitzOptions& opts,
                                    const core::TolerancePolicy& tol, SingularityInfo* info) {
  LuFactors fresh;
  if (!FactorizeMarkowitz(b, opts, tol.pivot, fresh, info)) return false;
  lu_ = std::move(fresh);
  etas_.clear();
  growth_ = 1.0;
  tol_ = tol;

  current_cols_.assign(b.cols, {});
  for (int c = 0; c < b.cols; ++c) {
    for (int p = b.col_ptr[c]; p < b.col_ptr[c + 1]; ++p) {
      current_cols_[c].emplace_back(b.row_idx[p], b.values[p]);
    }
  }
  return true;
}

std::vector<double> BasisFactorization::Ftran(const std::vector<std::pair<int, double>>& b_sparse) const {
  std::vector<double> x = la::Ftran(lu_, b_sparse);
  // Apply the eta chain in the order it was created: B^{-1} = E_r^{-1} ...
  // E_1^{-1} B_base^{-1}, so B^{-1} b = E_r^{-1}(...(E_1^{-1}(B_base^{-1} b))...).
  for (const auto& eta : etas_) {
    double yt = x[eta.slot] / eta.pivot;
    for (const auto& [i, a] : eta.off_diag) x[i] -= a * yt;
    x[eta.slot] = yt;
  }
  return x;
}

std::vector<double> BasisFactorization::Btran(const std::vector<std::pair<int, double>>& c_sparse) const {
  std::vector<double> y(lu_.m, 0.0);
  for (const auto& [orig_col, val] : c_sparse) y[orig_col] += val;

  // B^{-T} = B_base^{-T} E_1^{-T} ... E_r^{-T}, so apply the eta chain in
  // *reverse* creation order before the base BTRAN.
  for (auto it = etas_.rbegin(); it != etas_.rend(); ++it) {
    const EtaUpdate& eta = *it;
    double s = y[eta.slot];
    for (const auto& [i, a] : eta.off_diag) s -= a * y[i];
    y[eta.slot] = s / eta.pivot;
  }

  std::vector<std::pair<int, double>> y_sparse;
  y_sparse.reserve(y.size());
  for (int i = 0; i < static_cast<int>(y.size()); ++i) {
    if (y[i] != 0.0) y_sparse.emplace_back(i, y[i]);
  }
  return la::Btran(lu_, y_sparse);
}

bool BasisFactorization::Update(int t, const std::vector<std::pair<int, double>>& a_new_sparse,
                                 const core::TolerancePolicy& tol) {
  std::vector<double> alpha = Ftran(a_new_sparse);  // against the CURRENT basis, before this update
  double pivot = alpha[t];
  if (std::abs(pivot) < tol.pivot) return false;

  EtaUpdate eta;
  eta.slot = t;
  eta.pivot = pivot;
  double max_abs = std::abs(pivot);
  for (int i = 0; i < static_cast<int>(alpha.size()); ++i) {
    if (i == t || alpha[i] == 0.0) continue;
    eta.off_diag.emplace_back(i, alpha[i]);
    max_abs = std::max(max_abs, std::abs(alpha[i]));
  }
  growth_ = std::max(growth_, max_abs / std::abs(pivot));
  etas_.push_back(std::move(eta));

  current_cols_[t].assign(a_new_sparse.begin(), a_new_sparse.end());
  return true;
}

std::vector<double> BasisFactorization::MultiplyDense(const std::vector<double>& x) const {
  std::vector<double> y(lu_.m, 0.0);
  for (int c = 0; c < static_cast<int>(current_cols_.size()); ++c) {
    double xc = x[c];
    if (xc == 0.0) continue;
    for (const auto& [r, v] : current_cols_[c]) y[r] += v * xc;
  }
  return y;
}

std::vector<double> BasisFactorization::FtranRefined(const std::vector<std::pair<int, double>>& b_sparse,
                                                       int refinement_passes) const {
  std::vector<double> x = Ftran(b_sparse);
  std::vector<double> b_dense(lu_.m, 0.0);
  for (const auto& [i, v] : b_sparse) b_dense[i] += v;

  for (int pass = 0; pass < refinement_passes; ++pass) {
    std::vector<double> bx = MultiplyDense(x);
    std::vector<std::pair<int, double>> r_sparse;
    for (int i = 0; i < lu_.m; ++i) {
      double ri = b_dense[i] - bx[i];
      if (ri != 0.0) r_sparse.emplace_back(i, ri);
    }
    if (r_sparse.empty()) break;
    std::vector<double> dx = Ftran(r_sparse);
    for (int i = 0; i < lu_.m; ++i) x[i] += dx[i];
  }
  return x;
}

bool BasisFactorization::ShouldRefactorize() const {
  // How many eta updates to accumulate before rebuilding the factors from
  // scratch. This is a real trade-off, and it was measured rather than
  // guessed: a fresh factorization is expensive and its cost grows with m,
  // while every eta appended makes each subsequent FTRAN/BTRAN walk one
  // more link of the chain. Timing pilot87 (m ~ 2030, 22322 iterations)
  // end to end at several caps:
  //
  //     cap    total    refactorizations   share of time in refactorize
  //     100    95.7s          196                    77%
  //     200    47.5s           66                    59%
  //     400    41.1s           33                    34%
  //     800    did not finish inside 10 minutes
  //
  // A flat 100 spent three quarters of the solve rebuilding factors; 800
  // went off a cliff the other way, with the eta chain making every solve
  // against the basis ruinous. Scaling with m rather than picking one flat
  // number is the point: the bigger the basis, the more updates it takes
  // to be worth paying for a rebuild. The floor keeps small problems from
  // refactorizing constantly, and the ceiling keeps the chain bounded on
  // very large ones — 800 is proof that "just refactorize less" stops
  // being true well before the chain length stops growing.
  //
  // Numerical safety does not rest on this number: the growth monitor
  // below still forces a refactorization whenever the eta chain actually
  // becomes unstable, independent of how long it is.
  int cap = lu_.m / 10;
  if (cap < 50) cap = 50;
  if (cap > 250) cap = 250;
  if (static_cast<int>(etas_.size()) >= cap) return true;
  if (growth_ >= tol_.growth_refactor) return true;
  return false;
}

}  // namespace inferno::la
