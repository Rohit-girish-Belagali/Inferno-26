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

const char* NodeOutcomeName(NodeOutcome outcome) {
  switch (outcome) {
    case NodeOutcome::kRootRelaxation: return "root";
    case NodeOutcome::kBranched: return "branched";
    case NodeOutcome::kIntegerFeasible: return "integer";
    case NodeOutcome::kInfeasible: return "infeasible";
    case NodeOutcome::kDominated: return "dominated";
    case NodeOutcome::kGapCut: return "gap-cut";
    case NodeOutcome::kRelaxationFailed: return "unresolved";
  }
  return "unknown";
}

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
      core::Solution reduced = simplex::SolveRevised(pre.reduced, opts.lp_iteration_limit, tol);
      if (reduced.status != core::SolveStatus::kOptimal) {
        s.status = reduced.status;
        return reduced.status == core::SolveStatus::kInfeasible ? NodeResult::kInfeasible
                                                                 : NodeResult::kFailed;
      }
      s = presolve::Postsolve(lp, pre, reduced);
    } else {
      s = simplex::SolveRevised(lp, opts.lp_iteration_limit, tol);
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
  // Tries one specific rounding of a relaxation and keeps it if feasible
  // and better. `mode` 0 rounds to nearest, 1 rounds DOWN.
  //
  // Rounding down matters more than it looks. These models are dominated
  // by packing rows (sum of weights <= capacity), and rounding a
  // fractional solution UP almost always violates them -- which is why
  // nearest-rounding alone found NO incumbent at all on the 1000-binary
  // scale instance. Rounding down sacrifices objective value but lands
  // inside the feasible region, and any incumbent is worth far more than
  // none: without one, bound-based pruning cannot fire and the search
  // explores blindly.
  auto try_rounding_mode = [&](const std::vector<double>& x, int mode, double& incumbent_obj,
                               std::vector<double>& incumbent_x) {
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
        double target = (mode == 0) ? std::round(r[j]) : std::floor(r[j] + tol.feasibility);
        r[j] = std::min(std::max(target, lo), hi);
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

  auto try_rounding = [&](const std::vector<double>& x, double& incumbent_obj,
                          std::vector<double>& incumbent_x) {
    if (!opts.rounding_heuristic) return;
    try_rounding_mode(x, 0, incumbent_obj, incumbent_x);  // nearest
    try_rounding_mode(x, 1, incumbent_obj, incumbent_x);  // down: feasible for packing rows
  };

  // --- Observer plumbing. Everything reported here is read out of the
  // search state at the moment a node is disposed of; nothing is
  // interpolated, smoothed or predicted. When no callback is installed
  // none of this runs, so the search pays nothing for the capability.
  int events_emitted = 0;
  double incumbent_for_events = kInfinity;
  double bound_for_events = -kInfinity;
  const std::vector<Node>* stack_for_events = nullptr;
  auto emit = [&](int node_index, int depth, int branch_var, double branch_value, double node_bound,
                  NodeOutcome outcome) {
    if (!opts.node_callback) return;
    // The event budget bounds volume on a large tree. An incumbent
    // improvement is always reported regardless -- it is the one event a
    // watcher cannot afford to miss, and there are at most as many of them
    // as there are improvements.
    const bool always = (outcome == NodeOutcome::kIntegerFeasible);
    if (!always && events_emitted >= opts.event_limit) return;
    ++events_emitted;
    NodeEvent ev;
    ev.node_index = node_index;
    ev.depth = depth;
    ev.branch_var = branch_var;
    ev.branch_value = branch_value;
    ev.node_bound = node_bound;
    ev.incumbent = incumbent_for_events;
    // A valid global lower bound right now: the incumbent, every subtree
    // still open, and every subtree abandoned unresolved. Never reported
    // tighter than that minimum, so a watcher reading it mid-solve is
    // reading a real bound rather than an optimistic one.
    double lb = std::min(incumbent_for_events, bound_for_events);
    if (stack_for_events != nullptr) {
      for (const auto& nd : *stack_for_events) lb = std::min(lb, nd.parent_bound);
    }
    ev.best_bound = lb;
    ev.outcome = outcome;
    ev.elapsed_seconds = Now() - t_start;
    ev.open_nodes = stack_for_events ? static_cast<int>(stack_for_events->size()) : 0;
    opts.node_callback(ev);
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
    // Same reasoning as the no-incumbent path below: with no usable root
    // relaxation there is no objective, no bound and no gap, and leaving
    // the zero defaults in place would print like a perfectly solved
    // problem. This path is reached when the root LP itself exceeds the
    // iteration budget, which is exactly what happens on the largest
    // scale instances.
    out.objective_value = kInfinity;
    out.best_bound = -kInfinity;
    out.gap = kInfinity;
    out.proved_optimal = false;
    emit(0, 0, -1, 0.0, 0.0,
         root_res == NodeResult::kInfeasible ? NodeOutcome::kInfeasible
                                              : NodeOutcome::kRelaxationFailed);
    return out;
  }

  double incumbent_obj = kInfinity;
  std::vector<double> incumbent_x;
  try_rounding(root_sol.x, incumbent_obj, incumbent_x);
  incumbent_for_events = incumbent_obj;
  bound_for_events = root_sol.objective_value;
  emit(0, 0, -1, 0.0, root_sol.objective_value, NodeOutcome::kRootRelaxation);

  std::vector<Node> stack;
  root.parent_bound = root_sol.objective_value;
  stack.push_back(root);

  stack_for_events = &stack;
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
    if (node.parent_bound >= incumbent_obj - 1e-9) {
      emit(nodes, node.depth, -1, 0.0, node.parent_bound, NodeOutcome::kDominated);
      continue;
    }

    // Optional gap-based cut. When enabled this is a real trade: it stops
    // early, and it is recorded as unresolved so the reported bound and
    // proved_optimal both tell the truth about it.
    if (opts.gap_tolerance > 0.0 && std::isfinite(incumbent_obj)) {
      double allowed = opts.gap_tolerance * (std::abs(incumbent_obj) + 1e-10);
      if (node.parent_bound >= incumbent_obj - allowed) {
        mark_unresolved(node.parent_bound);
        bound_for_events = std::min(bound_for_events, best_unresolved);
        emit(nodes, node.depth, -1, 0.0, node.parent_bound, NodeOutcome::kGapCut);
        continue;
      }
    }

    core::Solution s;
    ++nodes;
    NodeResult res = solve_relaxation(node, s);
    if (res == NodeResult::kInfeasible) {
      emit(nodes, node.depth, -1, 0.0, node.parent_bound, NodeOutcome::kInfeasible);
      continue;  // proven empty: resolved
    }
    if (res == NodeResult::kFailed) {
      // We could not evaluate this subtree. Discarding it silently would
      // be the classic way to report a confident wrong bound.
      mark_unresolved(node.parent_bound);
      bound_for_events = std::min(bound_for_events, best_unresolved);
      emit(nodes, node.depth, -1, 0.0, node.parent_bound, NodeOutcome::kRelaxationFailed);
      if (opts.verbose) printf("  [node %d] relaxation failed; subtree left unresolved\n", nodes);
      continue;
    }

    if (s.objective_value >= incumbent_obj - 1e-9) {
      emit(nodes, node.depth, -1, 0.0, s.objective_value, NodeOutcome::kDominated);
      continue;  // proven dominated
    }

    int branch_var = most_fractional(s.x);
    if (branch_var == -1) {
      double v = s.objective_value;
      if (v < incumbent_obj) {
        incumbent_obj = v;
        incumbent_x = s.x;
        if (opts.verbose) printf("  [node %d] incumbent %.10g\n", nodes, incumbent_obj);
      }
      incumbent_for_events = incumbent_obj;
      emit(nodes, node.depth, -1, 0.0, s.objective_value, NodeOutcome::kIntegerFeasible);
      continue;
    }

    try_rounding(s.x, incumbent_obj, incumbent_x);
    incumbent_for_events = incumbent_obj;

    double v = s.x[branch_var];
    Node down = node, up = node;
    down.parent_bound = up.parent_bound = s.objective_value;
    down.depth = up.depth = node.depth + 1;
    down.upper.emplace_back(branch_var, std::floor(v));
    up.lower.emplace_back(branch_var, std::ceil(v));
    stack.push_back(up);
    stack.push_back(down);
    emit(nodes, node.depth, branch_var, v, s.objective_value, NodeOutcome::kBranched);
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
    // With no incumbent there is no objective and no gap. Leaving the
    // defaults in place reports objective 0 and gap 0, which reads exactly
    // like a perfectly solved problem -- the most misleading possible
    // output for a run that found nothing. Say "unknown" instead.
    out.objective_value = kInfinity;
    out.gap = kInfinity;
    out.proved_optimal = false;
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
