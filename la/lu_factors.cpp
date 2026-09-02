#include "la/lu_factors.hpp"

namespace inferno::la {

void BuildTransposeViews(LuFactors& lu) {
  int m = lu.m;
  lu.l_row.assign(m, {});
  lu.u_col.assign(m, {});

  for (int k = 0; k < m; ++k) {
    for (const auto& [kp, mult] : lu.l_col[k]) {
      // l_col[k] holds (k', mult) with k' > k, i.e. L[k'][k] = mult.
      // Its transpose entry lives at l_row[k'] as (k, mult) = L^T[k][k'']...
      // concretely: row view of L at row k' gets a (k, mult) entry (k < k').
      lu.l_row[kp].emplace_back(k, mult);
    }
    for (const auto& [kp, val] : lu.u_row[k]) {
      // u_row[k] holds (k', val) with k' > k, i.e. U[k][k'] = val.
      // Column view of U at column k' gets a (k, val) entry (k < k').
      lu.u_col[kp].emplace_back(k, val);
    }
  }
}

}  // namespace inferno::la
