#include "la/basis_factorization.hpp"

#include <algorithm>
#include <cmath>

#include "la/lu_solve.hpp"

namespace inferno::la {

bool BasisFactorization::Factorize(const core::CscMatrix& b, const MarkowitzOptions& opts,
                                    const core::TolerancePolicy& tol) {
  LuFactors fresh;
  if (!FactorizeMarkowitz(b, opts, tol.pivot, fresh)) return false;
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
  constexpr int kMaxUpdates = 100;
  if (static_cast<int>(etas_.size()) >= kMaxUpdates) return true;
  if (growth_ >= tol_.growth_refactor) return true;
  return false;
}

}  // namespace inferno::la
