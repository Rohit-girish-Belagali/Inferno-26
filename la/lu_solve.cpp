#include "la/lu_solve.hpp"

#include <algorithm>

namespace inferno::la {

namespace {

// DFS reachability over an adjacency list where adj[k] lists the indices
// that depend on k (i.e. must be visited/updated once k is finalized).
// Returns the reach set unsorted; callers sort ascending or descending
// depending on which triangular direction they're processing.
std::vector<int> ComputeReach(const std::vector<std::vector<std::pair<int, double>>>& adj,
                               const std::vector<int>& seeds, int m) {
  std::vector<char> visited(m, 0);
  std::vector<int> reach;
  std::vector<int> stack;
  for (int s : seeds) {
    if (s < 0 || s >= m || visited[s]) continue;
    visited[s] = 1;
    stack.push_back(s);
  }
  while (!stack.empty()) {
    int k = stack.back();
    stack.pop_back();
    reach.push_back(k);
    for (const auto& [kp, val] : adj[k]) {
      if (!visited[kp]) {
        visited[kp] = 1;
        stack.push_back(kp);
      }
    }
  }
  return reach;
}

}  // namespace

std::vector<double> Ftran(const LuFactors& lu, const std::vector<std::pair<int, double>>& b_sparse) {
  int m = lu.m;
  std::vector<double> y(m, 0.0);
  std::vector<int> seeds;
  seeds.reserve(b_sparse.size());
  for (const auto& [orig_row, val] : b_sparse) {
    int k = lu.row_perm_inv[orig_row];
    y[k] += val;
    seeds.push_back(k);
  }

  // Phase 1: L y = P b (forward, ascending k; right-looking via l_col,
  // whose entries always point to strictly larger k).
  std::vector<int> reach1 = ComputeReach(lu.l_col, seeds, m);
  std::sort(reach1.begin(), reach1.end());
  for (int k : reach1) {
    double yk = y[k];
    if (yk == 0.0) continue;
    for (const auto& [kp, mult] : lu.l_col[k]) y[kp] -= mult * yk;
  }

  // Phase 2: U z = y (backward, descending k; right-looking via u_col,
  // whose entries always point to strictly smaller k), result left in y.
  std::vector<int> seeds2;
  for (int k : reach1) {
    if (y[k] != 0.0) seeds2.push_back(k);
  }
  std::vector<int> reach2 = ComputeReach(lu.u_col, seeds2, m);
  std::sort(reach2.rbegin(), reach2.rend());
  for (int k : reach2) {
    double zk = y[k] / lu.diag[k];
    y[k] = zk;
    if (zk == 0.0) continue;
    for (const auto& [kp, val] : lu.u_col[k]) y[kp] -= val * zk;
  }

  std::vector<double> x(m, 0.0);
  for (int k = 0; k < m; ++k) x[lu.col_perm[k]] = y[k];
  return x;
}

std::vector<double> Btran(const LuFactors& lu, const std::vector<std::pair<int, double>>& c_sparse) {
  int m = lu.m;
  std::vector<double> v(m, 0.0);
  std::vector<int> seeds;
  seeds.reserve(c_sparse.size());
  for (const auto& [orig_col, val] : c_sparse) {
    int k = lu.col_perm_inv[orig_col];
    v[k] += val;
    seeds.push_back(k);
  }

  // Phase 1: U^T w = c (forward, ascending k; right-looking via u_row,
  // whose entries always point to strictly larger k — U^T[k'][k] = U[k][k']
  // for k' > k is exactly u_row[k]).
  std::vector<int> reach1 = ComputeReach(lu.u_row, seeds, m);
  std::sort(reach1.begin(), reach1.end());
  for (int k : reach1) {
    double wk = v[k] / lu.diag[k];
    v[k] = wk;
    if (wk == 0.0) continue;
    for (const auto& [kp, val] : lu.u_row[k]) v[kp] -= val * wk;
  }

  // Phase 2: L^T u = w (backward, descending k; right-looking via l_row,
  // whose entries always point to strictly smaller k; unit diagonal).
  std::vector<int> seeds2;
  for (int k : reach1) {
    if (v[k] != 0.0) seeds2.push_back(k);
  }
  std::vector<int> reach2 = ComputeReach(lu.l_row, seeds2, m);
  std::sort(reach2.rbegin(), reach2.rend());
  for (int k : reach2) {
    double uk = v[k];
    if (uk == 0.0) continue;
    for (const auto& [kp, mult] : lu.l_row[k]) v[kp] -= mult * uk;
  }

  std::vector<double> y(m, 0.0);
  for (int k = 0; k < m; ++k) y[lu.row_perm[k]] = v[k];
  return y;
}

}  // namespace inferno::la
