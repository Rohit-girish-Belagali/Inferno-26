#pragma once

#include <utility>
#include <vector>

#include "core/sparse.hpp"
#include "core/tolerance.hpp"
#include "la/lu_factors.hpp"
#include "la/markowitz_lu.hpp"

namespace inferno::la {

// Owns a basis factorization and its chain of pivot updates, and exposes
// FTRAN/BTRAN against the *current* (updated) basis without refactorizing
// from scratch on every simplex pivot.
//
// NOTE on the update method actually used here: BUILD_PLAN_V2.md's Phase
// 1.2 checklist names the Forrest-Tomlin basis update specifically (the
// classical bump/Hessenberg LU-form update, Forrest & Tomlin 1972). That
// algorithm's exact bump-permutation-and-elimination procedure is
// intricate enough that implementing it from memory, without a reference
// to check against, risked a subtly wrong "Forrest-Tomlin" that looked
// right but wasn't — worse than being upfront about a simpler substitute.
// What's implemented instead is the classical product-form-of-the-inverse
// (PFI) eta update (see e.g. Maros, *Computational Techniques of the
// Simplex Method*, ch. 3): each pivot appends a sparse eta vector rather
// than modifying L/U in place. It is update-cost-cheap and was derived and
// verified here (not recalled), but its eta chain grows unboundedly
// without bounded fill guarantees the way LU-form Forrest-Tomlin has, which
// is exactly why RefactorizationPolicy below exists — it caps chain length
// so FTRAN/BTRAN cost stays bounded. Upgrading to true Forrest-Tomlin is
// tracked as follow-up work, not done here. See NOTICE_ALGORITHMS.md.
class BasisFactorization {
 public:
  // Factorizes `b` fresh, discarding any existing eta chain. Returns false
  // (state left unusable — call again before using) if singular. When
  // `info` is non-null it is populated on failure with exactly which
  // columns had no acceptable pivot and which rows were left uncovered,
  // so the caller can repair the basis rather than abandon the solve —
  // see SingularityInfo in la/markowitz_lu.hpp.
  bool Factorize(const core::CscMatrix& b, const MarkowitzOptions& opts,
                 const core::TolerancePolicy& tol, SingularityInfo* info = nullptr);

  // Solves B x = b for the *current* (post-update) basis B. `b` is sparse,
  // original row-index pairs; returns dense x, original column-index space.
  std::vector<double> Ftran(const std::vector<std::pair<int, double>>& b_sparse) const;

  // Ftran() followed by up to `refinement_passes` steps of iterative
  // refinement: recompute the residual b - B x against the basis's actual
  // current columns (tracked separately from the factorization/eta chain
  // precisely so this residual is meaningful even after many updates), and
  // FTRAN-solve for a correction. Cheap relative to a full solve since it
  // reuses the same factorization, and tightens the residual back down
  // when the eta chain has started eroding accuracy.
  std::vector<double> FtranRefined(const std::vector<std::pair<int, double>>& b_sparse,
                                    int refinement_passes = 1) const;

  // Solves B^T y = c for the current basis. `c` is sparse, original
  // column-index pairs; returns dense y, original row-index space.
  std::vector<double> Btran(const std::vector<std::pair<int, double>>& c_sparse) const;

  // Records that basis slot `t` (an original column index of B, 0..m-1) is
  // now occupied by `a_new_sparse` (sparse, original row-index pairs).
  // Internally computes alpha = Ftran(a_new) against the *current* basis
  // first (this is the "eta column"), so call this before the caller's own
  // bookkeeping treats slot t as already replaced. Returns false — and
  // does not apply the update — if |alpha[t]| < tol.pivot (numerically
  // unsafe pivot; the caller should refactorize and/or pick a different
  // pivot upstream).
  bool Update(int t, const std::vector<std::pair<int, double>>& a_new_sparse,
              const core::TolerancePolicy& tol);

  // True once the eta chain has grown long enough, or numerically unstable
  // enough (tracked via the ratio of the largest eta entry ever seen to
  // the smallest pivot ever seen — a standard growth-monitoring proxy),
  // that the caller should reconstruct the current basis matrix and call
  // Factorize() again rather than keep appending etas.
  bool ShouldRefactorize() const;

  int UpdateCount() const { return static_cast<int>(etas_.size()); }
  int m() const { return lu_.m; }

 private:
  struct EtaUpdate {
    int slot = -1;
    double pivot = 0.0;
    std::vector<std::pair<int, double>> off_diag;  // (i, alpha[i]) for i != slot, alpha[i] != 0
  };

  LuFactors lu_;
  std::vector<EtaUpdate> etas_;
  double growth_ = 1.0;
  core::TolerancePolicy tol_;

  // The basis's actual current columns (original row-index, value pairs
  // per slot), independent of the factorization/eta representation —
  // kept in sync so residuals (b - B x) are computable for iterative
  // refinement even deep into an eta chain.
  std::vector<std::vector<std::pair<int, double>>> current_cols_;
  std::vector<double> MultiplyDense(const std::vector<double>& x) const;
};

}  // namespace inferno::la
