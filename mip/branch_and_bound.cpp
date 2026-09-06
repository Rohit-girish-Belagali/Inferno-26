#include "mip/branch_and_bound.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "presolve/presolve.hpp"
#include "simplex/revised_simplex.hpp"

namespace inferno::mip {

namespace {

using core::kInfinity;

double Now() {
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// A node is just the bound changes that distinguish it from the root.
// Storing the deltas rather than a whole copy of the problem keeps memory
// proportional to depth instead of to the number of nodes, which is what
// makes a deep tree affordable.
struct Node {
  std::vector<std::pair<int, double>> lower;  // column -> tightened lower bound
  std::vector<std::pair<int, double>> upper;  // column -> tightened upper bound
  double parent_bound = -kInfinity;
  int depth = 0;
};

double Frac(double v) { return v - std::floor(v); }

}  // namespace

core::MipSolution SolveMip(const core::MipProblem& problem, const BranchAndBoundOptions& opts,
                            const core::TolerancePolicy& tol) {
  core::MipSolution out;
  const int n = problem.lp.num_cols;
  const double t_start = Now();

  auto solve_relaxation = [&](const Node& node, core::Solution& s) {
    core::LpProblem lp = problem.lp;
    for (const auto& [j, v] : node.lower) lp.col_lo[j] = std::max(lp.col_lo[j], v);
    for (const auto& [j, v] : node.upper) lp.col_hi[j] = std::min(lp.col_hi[j], v);
    // An empty box means the branch is infeasible by construction; say so
    // without troubling the simplex.
    for (int j = 0; j < n; ++j) {
      if (lp.col_lo[j] > lp.col_hi[j] + tol.feasibility) return false;
    }
    if (opts.use_presolve) {
      auto pre = presolve::Presolve(lp);
      if (pre.infeasible) return false;
      core::Solution reduced = simplex::SolveRevised(pre.reduced);
      if (reduced.status != core::SolveStatus::kOptimal) {
        s.status = reduced.status;
        return false;
      }
      s = presolve::Postsolve(lp, pre, reduced);
    } else {
      s = simplex::SolveRevised(lp);
    }
    return s.status == core::SolveStatus::kOptimal;
  };

  // Index of the "most fractional" integer variable, or -1 if the point is
  // already integral. Most-fractional is the simplest sound rule: it
  // branches where the relaxation is least decided, which tends to move
  // the bound on both children rather than only one.
  auto most_fractional = [&](const std::vector<double>& x) {
    int best = -1;
    double best_score = opts.integrality_tolerance;
    for (int j = 0; j < n; ++j) {
      if (!problem.is_integer[j]) continue;
      double f = Frac(x[j]);
      double score = std::min(f, 1.0 - f);  // distance to the nearest integer
      if (score > best_score) {
        best_score = score;
        best = j;
      }
    }
    return best;
  };

  auto objective_of = [&](const std::vector<double>& x) {
    double v = problem.lp.obj_offset;
    for (int j = 0; j < n; ++j) v += problem.lp.obj[j] * x[j];
    return v;
  };

  // Rounds a relaxation to the nearest integer point and keeps it if it is
  // actually feasible. Cheap, and an incumbent found early is what lets
  // bound pruning discard subtrees instead of exploring them.
  auto try_rounding = [&](const std::vector<double>& x, double& incumbent_obj,
                          std::vector<double>& incumbent_x) {
    if (!opts.rounding_heuristic) return;
    std::vector<double> r = x;
    for (int j = 0; j < n; ++j) {
      if (problem.is_integer[j]) r[j] = std::round(r[j]);
      r[j] = std::min(std::max(r[j], problem.lp.col_lo[j]), problem.lp.col_hi[j]);
      if (problem.is_integer[j]) r[j] = std::round(r[j]);
    }
    std::vector<double> ax(problem.lp.num_rows, 0.0);
    problem.lp.a.MultiplyAdd(r, ax);
    for (int i = 0; i < problem.lp.num_rows; ++i) {
      if (std::isfinite(problem.lp.row_lo[i]) && ax[i] < problem.lp.row_lo[i] - tol.feasibility) return;
      if (std::isfinite(problem.lp.row_hi[i]) && ax[i] > problem.lp.row_hi[i] + tol.feasibility) return;
    }
    double v = objective_of(r);
    if (v < incumbent_obj) {
      incumbent_obj = v;
      incumbent_x = r;
    }
  };

  // --- Root relaxation. Its objective is the initial global bound. ---
  Node root;
  core::Solution root_sol;
  if (!solve_relaxation(root, root_sol)) {
    out.status = (root_sol.status == core::SolveStatus::kUnbounded)
                     ? core::SolveStatus::kUnbounded
                     : core::SolveStatus::kInfeasible;
    return out;
  }

  double incumbent_obj = kInfinity;
  std::vector<double> incumbent_x;
  try_rounding(root_sol.x, incumbent_obj, incumbent_x);

  std::vector<Node> stack;
  root.parent_bound = root_sol.objective_value;
  stack.push_back(root);

  double global_bound = root_sol.objective_value;
  int nodes = 0;

  while (!stack.empty()) {
    if (nodes >= opts.node_limit || Now() - t_start > opts.time_limit_seconds) break;

    // Depth-first: take the most recently created node. Depth-first drives
    // toward a feasible integer solution quickly, and an incumbent is what
    // makes pruning effective; best-bound alone can explore enormously
    // before finding its first one.
    Node node = stack.back();
    stack.pop_back();

    // A node whose parent bound is already no better than the incumbent
    // cannot contain a strictly better solution. Strictly worse, by more
    // than a tolerance — never on equality, so a tie cannot silently lose
    // the optimum.
    if (node.parent_bound >= incumbent_obj - opts.gap_tolerance * std::abs(incumbent_obj) - 1e-9) {
      continue;
    }

    core::Solution s;
    ++nodes;
    if (!solve_relaxation(node, s)) continue;  // infeasible or unsolvable: prune

    if (s.objective_value >= incumbent_obj - 1e-9) continue;  // bound says no improvement here

    int branch_var = most_fractional(s.x);
    if (branch_var == -1) {
      // Integral: a genuine incumbent.
      double v = s.objective_value;
      if (v < incumbent_obj) {
        incumbent_obj = v;
        incumbent_x = s.x;
        if (opts.verbose) {
          printf("  [node %d] incumbent %.6g\n", nodes, incumbent_obj);
        }
      }
      continue;
    }

    try_rounding(s.x, incumbent_obj, incumbent_x);

    // Branch: x_j <= floor(v) and x_j >= ceil(v). Between them these cover
    // every integer value of x_j and exclude only the fractional interval,
    // so nothing integral is lost.
    double v = s.x[branch_var];
    Node down = node, up = node;
    down.parent_bound = up.parent_bound = s.objective_value;
    down.depth = up.depth = node.depth + 1;
    down.upper.emplace_back(branch_var, std::floor(v));
    up.lower.emplace_back(branch_var, std::ceil(v));
    // Push the "up" branch first so the "down" branch is explored first —
    // rounding down more often lands in the feasible region for the
    // capacity-style constraints these models are full of.
    stack.push_back(up);
    stack.push_back(down);
  }

  // The proven bound is the best relaxation value still unexplored; with
  // the stack empty the tree is exhausted and the incumbent is optimal.
  if (!stack.empty()) {
    double best_open = kInfinity;
    for (const auto& nd : stack) best_open = std::min(best_open, nd.parent_bound);
    global_bound = std::max(global_bound, std::min(best_open, incumbent_obj));
  } else {
    global_bound = incumbent_obj;
  }

  out.nodes_explored = nodes;
  out.best_bound = global_bound;
  if (incumbent_x.empty()) {
    out.status = core::SolveStatus::kInfeasible;
    return out;
  }
  out.x = incumbent_x;
  out.objective_value = incumbent_obj;
  out.gap = std::abs(incumbent_obj - global_bound) / (std::abs(incumbent_obj) + 1e-10);
  out.proved_optimal = stack.empty();
  out.status = core::SolveStatus::kOptimal;
  return out;
}

}  // namespace inferno::mip
