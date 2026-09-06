#pragma once

#include "core/qp_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::qp {

// ADMM for convex QP (BUILD_PLAN_V2.md Phase 4Q), completing the LP/QP
// pair of the three solver types the problem statement names.
//
// Splitting. The constraint set enters through a copy z = Ax, so the
// problem becomes
//
//     min 1/2 x'Px + q'x + I_C(z)   s.t.   Ax = z,     C = [row_lo, row_hi]
//
// which ADMM alternates: an equality-constrained quadratic solve in x
// (linear, and the same matrix every iteration), a projection of z onto
// the box (a clamp), and a dual update. Splitting is what buys that — the
// quadratic and the inequalities never have to be handled at once.
//
// The x-step is the KKT system
//
//     [ P + sigma I     A'    ] [ x ]   [ sigma x^k - q          ]
//     [ A            -1/rho I ] [ nu ] = [ z^k - y^k / rho       ]
//
// which is quasi-definite for sigma, rho > 0 — so it is nonsingular even
// when P is only positive SEMI-definite, and it is FIXED across iterations.
// That is the whole efficiency argument for ADMM here: factorize once,
// then every iteration is two triangular solves plus vector work. The
// factorization is Phase 1.2's Markowitz LU, reused as the plan intends
// rather than a second linear algebra stack written beside it.
//
// Convergence is linear at best, as for any first-order method — this
// reaches modest accuracy quickly and high accuracy slowly. Termination is
// on the KKT residuals actually measured, and `qp/checker` recomputes them
// independently afterwards.
struct AdmmOptions {
  int max_iterations = 20000;
  double rho = 0.1;      // penalty on the splitting constraint
  double sigma = 1e-6;   // regularization keeping the KKT system definite
  double alpha = 1.6;    // over-relaxation; 1.0 is plain ADMM
  int check_every = 25;
  // Rho is rebalanced when the primal and dual residuals drift apart, which
  // matters more than the starting value does — a fixed rho suits almost no
  // real problem, and refactorizing to change it is affordable precisely
  // because the factorization is amortized over many iterations.
  bool adaptive_rho = true;
};

core::QpSolution SolveQpAdmm(const core::QpProblem& problem,
                              const AdmmOptions& opts = AdmmOptions{},
                              const core::TolerancePolicy& tol = core::DefaultTolerances());

}  // namespace inferno::qp
