#pragma once

#include <functional>

#include "core/mip_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::mip {

// Branch-and-bound for mixed-integer linear programming.
//
// The method in one line: solve the LP relaxation; if a variable that
// must be integral comes back fractional, split the problem into two
// sub-problems that exclude that fractional value, and recurse. The LP
// relaxation's objective is a valid BOUND on every integer solution below
// it, which is what makes the search finite in practice — a node whose
// bound is already no better than the best integer solution found so far
// cannot contain a better one, so the entire subtree is discarded without
// being looked at. Everything else here exists to make that pruning
// happen sooner.
//
// Correctness rests on one invariant, and it is worth stating plainly
// because violating it is how MILP solvers return confident wrong
// answers: **pruning must never discard a subtree that could contain a
// strictly better integer solution than the incumbent.** The bound
// comparison below is therefore deliberately conservative — a node is cut
// only when its relaxation is worse than the incumbent by more than a
// tolerance, never on equality, so a tie can never silently lose the
// optimum.
// What actually happened at one node of the real search. This exists so a
// user interface can show the search as it happens WITHOUT the interface
// inventing anything: every field below is read straight out of the
// branch-and-bound state at the moment the node was disposed of. There is
// no simulated progress anywhere in this project, and this type is the
// reason there does not need to be.
enum class NodeOutcome {
  kRootRelaxation = 0,   // the root LP; its objective is the initial bound
  kBranched = 1,         // fractional: split into two children
  kIntegerFeasible = 2,  // relaxation came back integral
  kInfeasible = 3,       // PROVEN empty subtree
  kDominated = 4,        // bound proves it cannot beat the incumbent
  kGapCut = 5,           // stopped by gap tolerance: a trade, not a proof
  kRelaxationFailed = 6, // could not evaluate: left UNRESOLVED
};

const char* NodeOutcomeName(NodeOutcome outcome);

struct NodeEvent {
  int node_index = 0;    // 0 is the root relaxation
  int depth = 0;
  int branch_var = -1;   // column branched on, or -1
  double branch_value = 0.0;  // the fractional value that forced the branch
  double node_bound = 0.0;    // this node's relaxation objective
  double incumbent = 0.0;     // best integer objective so far (+inf if none)
  double best_bound = 0.0;    // valid global lower bound at this moment
  NodeOutcome outcome = NodeOutcome::kRootRelaxation;
  double elapsed_seconds = 0.0;
  int open_nodes = 0;    // size of the search stack after this node
};

using NodeCallback = std::function<void(const NodeEvent&)>;

struct BranchAndBoundOptions {
  int node_limit = 200000;
  double time_limit_seconds = 60.0;
  // Stop once the proven gap is this small. 0 demands a proof of
  // optimality; industrial practice usually accepts a fraction of a
  // percent, and the gap actually achieved is always reported.
  double gap_tolerance = 1e-6;
  // How far from an integer a value may sit and still count as integral.
  double integrality_tolerance = 1e-6;
  bool use_presolve = true;
  // Try rounding the relaxation to an integer point at each node. Cheap,
  // and an early incumbent is what makes bound-based pruning bite.
  bool rounding_heuristic = true;
  bool verbose = false;
  // Iteration cap for the LP solved at each node, INCLUDING the root.
  // Without it the time limit is not actually a limit: the root relaxation
  // is solved before the search loop begins and the loop's clock check
  // only runs between nodes, so a single slow LP runs to completion no
  // matter what the caller asked for. Measured: a 5000-binary instance
  // blew past a 15-second limit by more than ten minutes inside the root
  // solve alone.
  //
  // A node whose LP hits this cap is treated as UNRESOLVED, not as
  // infeasible -- which the bound accounting already handles correctly, so
  // capping work can cost a proof but can never produce a wrong answer.
  int lp_iteration_limit = 200000;

  // Optional observer, invoked once per disposed node with the real search
  // state. Installing one costs a running minimum over the open stack (to
  // report a valid global bound at each step), which is why it is only
  // computed when somebody is actually listening.
  NodeCallback node_callback;
  // Bound on how many events are emitted, so a million-node search does not
  // drown its observer. Once the budget is spent only incumbent
  // improvements and the final node continue to be reported -- the events
  // stop being exhaustive, they never start being fabricated.
  int event_limit = 4000;
};

core::MipSolution SolveMip(const core::MipProblem& problem,
                            const BranchAndBoundOptions& opts = BranchAndBoundOptions{},
                            const core::TolerancePolicy& tol = core::DefaultTolerances());

}  // namespace inferno::mip
