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

  // A node's relaxation can end three ways, and conflating them is how a
  // branch-and-bound reports a bound that is mathematically wrong:
  //   kOptimal    - usable bound, branch or accept
  //   kInfeasible - PROVEN empty, safe to discard entirely
  //   kFailed     - numerical error, unbounded, or iteration limit. We do
  //                 NOT know the subtree is empty, so discarding it is
  //                 unsound: the optimum may be inside. It is recorded as
  //                 unresolved instead, which both weakens the reported
  //                 bound and forbids claiming proved optimality.
  enum class NodeResult { kOptimal, kInfeasible, kFailed };

  auto solve_relaxation = [&](const Node& node, core::Solution& s) -> NodeResult {
    core::LpProblem lp = problem.lp;
    for (const auto& [j, v] : node.lower) lp.col_lo[j] = std::max(lp.col_lo[j], v);
    for (const auto& [j, v] : node.upper) lp.col_hi[j] = std::min(lp.col_hi[j], v);
    // An empty box means the branch is infeasible by construction; say so
    // without troubling the simplex.
    for (int j = 0; j < n; ++j) {
      if (lp.col_lo[j] > lp.col_hi[j] + tol.feasibility) return NodeResult::kInfeasible;
    }
    if (opts.use_presolve) {
      auto pre = presolve::Presolve(lp);
      if (pre.infeasible) return NodeResult::kInfeasible;
      core::Solution reduced = simplex::SolveRevised(pre.reduced);
      if (reduced.status != core::SolveStatus::kOptimal) {
        s.status = reduced.status;
        return reduced.status == core::SolveStatus::kInfeasible ? NodeResult::kInfeasible
                                                                 : NodeResult::kFailed;
      }
      s = presolve::Postsolve(lp, pre, reduced);
    } else {
      s = simplex::SolveRevised(lp);
    }
    if (s.status == core::SolveStatus::kOptimal) return NodeResult::kOptimal;
    if (s.status == core::SolveStatus::kInfeasible) return NodeResult::kInfeasible;
    return NodeResult::kFailed;
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
      if (problem.is_integer[j]) {
        // Round to the nearest integer that is actually INSIDE the column
        // box. Rounding after a clamp is not enough: with a fractional
        // bound (say col_hi = 2.5) it rounds straight back out to 3, and
        // the row-feasibility check below never looks at column bounds, so
        // an out-of-box point could become the incumbent and be reported
        // as the answer.
        double lo = std::isfinite(problem.lp.col_lo[j]) ? std::ceil(problem.lp.col_lo[j] - tol.feasibility)
                                                         : -kInfinity;
        double hi = std::isfinite(problem.lp.col_hi[j]) ? std::floor(problem.lp.col_hi[j] + tol.feasibility)
                                                         : kInfinity;
        if (lo > hi) return;  // no integer fits in this column's box
        r[j] = std::min(std::max(std::round(r[j]), lo), hi);
      } else {
        r[j] = std::min(std::max(r[j], problem.lp.col_lo[j]), problem.lp.col_hi[j]);
      }
    }
    // Verify the column box explicitly rather than assuming the clamping
    // above got it right.
    for (int j = 0; j < n; ++j) {
      if (std::isfinite(problem.lp.col_lo[j]) && r[j] < problem.lp.col_lo[j] - tol.feasibility) return;
      if (std::isfinite(problem.lp.col_hi[j]) && r[j] > problem.lp.col_hi[j] + tol.feasibility) return;
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
  NodeResult root_res = solve_relaxation(root, root_sol);
  if (root_res != NodeResult::kOptimal) {
    // The root's three outcomes are genuinely different answers about the
    // MILP and must not be collapsed. An infeasible relaxation proves the
    // MILP infeasible (the integer set is a subset of it); an unbounded
    // relaxation is reported as such; anything else is a failure to
    // determine, which is neither.
    if (root_res == NodeResult::kInfeasible) {
      out.status = core::SolveStatus::kInfeasible;
    } else if (root_sol.status == core::SolveStatus::kUnbounded) {
      out.status = core::SolveStatus::kUnbounded;
    } else {
      out.status = core::SolveStatus::kNumericalError;
    }
    return out;
  }

  double incumbent_obj = kInfinity;
  std::vector<double> incumbent_x;
  try_rounding(root_sol.x, incumbent_obj, incumbent_x);

  std::vector<Node> stack;
  root.parent_bound = root_sol.objective_value;
  stack.push_back(root);

  const double root_bound = root_sol.objective_value;
  int nodes = 0;

  // Bound accounting, rebuilt so it is valid by construction rather than
  // reconstructed at the end.
  //
  // For a minimisation, the optimum is either the incumbent or lies in
  // some subtree we did not fully resolve. So a valid lower bound is
  //     LB = min(incumbent, min over UNRESOLVED subtrees of their bound)
  // and a subtree counts as resolved only when it was explored to
  // completion, PROVEN infeasible, or PROVEN unable to beat the incumbent.
  // Anything else -- a failed relaxation, a node left on the stack when a
  // limit fired, a node cut by the gap tolerance -- is unresolved and must
  // weaken the bound. `all_resolved` is what licenses claiming optimality;
  // it is false the moment any subtree is abandoned for a reason other
  // than proof.
  double best_unresolved = kInfinity;
  bool all_resolved = true;
  auto mark_unresolved = [&](double bound) {
    best_unresolved = std::min(best_unresolved, bound);
    all_resolved = false;
  };

  while (!stack.empty()) {
    if (nodes >= opts.node_limit || Now() - t_start > opts.time_limit_seconds) break;

    Node node = stack.back();
    stack.pop_back();

    // Bound dominance. A subtree whose bound is already >= the incumbent
    // cannot contain anything strictly better, so discarding it is a
    // PROOF, not an approximation, and the bound is untouched. The
    // comparison is strict-with-epsilon and deliberately does NOT use
    // gap_tolerance: cutting on the gap tolerance would discard subtrees
    // that might genuinely be better, which is fine as a speed/accuracy
    // trade but must not then be reported as proved optimality.
    if (node.parent_bound >= incumbent_obj - 1e-9) continue;

    // Optional gap-based cut. When enabled this is a real trade: it stops
    // early, and it is recorded as unresolved so the reported bound and
    // proved_optimal both tell the truth about it.
    if (opts.gap_tolerance > 0.0 && std::isfinite(incumbent_obj)) {
      double allowed = opts.gap_tolerance * (std::abs(incumbent_obj) + 1e-10);
      if (node.parent_bound >= incumbent_obj - allowed) {
        mark_unresolved(node.parent_bound);
        continue;
      }
    }

    core::Solution s;
    ++nodes;
    NodeResult res = solve_relaxation(node, s);
    if (res == NodeResult::kInfeasible) continue;  // proven empty: resolved
    if (res == NodeResult::kFailed) {
      // We could not evaluate this subtree. Discarding it silently would
      // be the classic way to report a confident wrong bound.
      mark_unresolved(node.parent_bound);
      if (opts.verbose) printf("  [node %d] relaxation failed; subtree left unresolved\n", nodes);
      continue;
    }

    if (s.objective_value >= incumbent_obj - 1e-9) continue;  // proven dominated

    int branch_var = most_fractional(s.x);
    if (branch_var == -1) {
      double v = s.objective_value;
      if (v < incumbent_obj) {
        incumbent_obj = v;
        incumbent_x = s.x;
        if (opts.verbose) printf("  [node %d] incumbent %.10g\n", nodes, incumbent_obj);
      }
      continue;
    }

    try_rounding(s.x, incumbent_obj, incumbent_x);

    double v = s.x[branch_var];
    Node down = node, up = node;
    down.parent_bound = up.parent_bound = s.objective_value;
    down.depth = up.depth = node.depth + 1;
    down.upper.emplace_back(branch_var, std::floor(v));
    up.lower.emplace_back(branch_var, std::ceil(v));
    stack.push_back(up);
    stack.push_back(down);
  }

  // Whatever is still on the stack when a limit fired is unresolved.
  for (const auto& nd : stack) mark_unresolved(nd.parent_bound);

  // LB = min(incumbent, best unresolved bound), never below the root
  // relaxation, which is itself always a valid bound.
  double global_bound = std::min(incumbent_obj, best_unresolved);
  if (!std::isfinite(global_bound)) global_bound = root_bound;
  global_bound = std::max(global_bound, root_bound);

  out.nodes_explored = nodes;
  out.best_bound = global_bound;
  if (incumbent_x.empty()) {
    // No integer point found. That is only a PROOF of infeasibility if the
    // whole tree was resolved; otherwise the honest answer is that the
    // search ran out of budget without finding one.
    out.status = all_resolved ? core::SolveStatus::kInfeasible : core::SolveStatus::kIterationLimit;
    return out;
  }
  out.x = incumbent_x;
  out.objective_value = incumbent_obj;
  out.gap = std::abs(incumbent_obj - global_bound) / (std::abs(incumbent_obj) + 1e-10);
  // Optimality is claimed only when every subtree was actually resolved.
  // stack.empty() alone is not enough: nodes dropped for a failed
  // relaxation or cut by the gap tolerance leave the stack empty while
  // proving nothing.
  out.proved_optimal = all_resolved;
  out.status = core::SolveStatus::kOptimal;
  return out;
}

}  // namespace inferno::mip
