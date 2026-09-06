#pragma once

#include "core/lp_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::firstorder {

// PDLP — a primal-dual first-order LP solver (BUILD_PLAN_V2.md Phase 3).
//
// Why this exists alongside a working simplex: simplex walks vertex to
// vertex along a path and cannot be parallelized. PDLP is a fixed-point
// iteration built from two sparse matrix-vector products and vector
// arithmetic, so it is bandwidth-bound and maps onto a GPU. This is the
// CPU reference implementation the plan asks for first; the CUDA port is
// blocked on hardware (the dev machine is an Apple M4 with no NVIDIA GPU
// and no nvcc), which is recorded in the plan as a known blocker rather
// than discovered late.
//
// The problem, as this project stores it:
//
//     min  c'x   s.t.   row_lo <= Ax <= row_hi,   col_lo <= x <= col_hi
//
// Written with indicator functions, min_x c'x + I_X(x) + I_C(Ax), where
// X is the column box and C the row box. Its saddle-point form is
//
//     min_x max_y  c'x + I_X(x) + y'(Ax) - h*(y),      h = I_C
//
// and PDHG (Chambolle-Pock) alternates a proximal ascent step in y with a
// proximal descent step in x:
//
//     y+ = prox_{sigma h*}( y + sigma A xbar )
//     x+ = proj_X( x - tau (c + A' y+) )
//     xbar+ = 2 x+ - x
//
// The dual prox is the one piece worth spelling out, because h* is the
// support function of the row box and has no convenient closed form on its
// own. Moreau decomposition turns it into a projection onto C instead:
//
//     prox_{sigma h*}(v) = v - sigma * proj_C(v / sigma)
//
// which is a clamp, so the whole iteration is two SpMVs plus vector work.
//
// Convergence requires tau * sigma * ||A||^2 < 1; ||A||_2 is estimated by
// power iteration up front. The ratio between tau and sigma is the
// "primal weight" and is rebalanced during the solve, since a single fixed
// ratio suits almost no real instance.
//
// HONEST SCOPE. A first-order method converges linearly at best and is
// slow to high accuracy — it buys scalability and parallelism, not the
// crisp 1e-9 vertex solutions simplex produces. This returns kOptimal only
// when the relative KKT residuals it measures are actually below
// `tol.checker_residual`, and the independent checker still has the final
// say exactly as it does for every other solve path here. What is NOT
// implemented yet, and is not pretended otherwise: the CUDA port, and
// crossover to a vertex solution.
struct PdlpOptions {
  int max_iterations = 100000;
  // How often to evaluate the KKT residuals and consider a restart. The
  // check costs two SpMVs, so doing it every iteration would roughly halve
  // the throughput of an iteration that is otherwise just SpMV + vector
  // work.
  int check_every = 64;
  bool use_restarts = true;
  bool use_primal_weight_balancing = true;
};

struct PdlpResult {
  core::Solution solution;
  int iterations = 0;
  double relative_primal_residual = 0.0;
  double relative_dual_residual = 0.0;
  double relative_duality_gap = 0.0;
};

PdlpResult SolvePdlp(const core::LpProblem& problem, const PdlpOptions& opts = PdlpOptions{},
                      const core::TolerancePolicy& tol = core::DefaultTolerances());

}  // namespace inferno::firstorder
