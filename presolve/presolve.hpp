#pragma once

#include <utility>
#include <vector>

#include "core/lp_problem.hpp"
#include "core/tolerance.hpp"

namespace inferno::presolve {

// Phase 2.2 — shrink the problem before solving, then exactly reverse
// every reduction to recover a solution to the ORIGINAL problem.
//
// Six reduction kinds, applied to a fixpoint (fixing one variable can make
// another variable's row bounds collapse it to fixed too; removing a row
// can make a column singleton in another row, etc.):
//
//   1. Fixed variable (col_lo == col_hi): substitute the known value into
//      every row it appears in and into the objective, then remove the
//      column. No dual recovery needed (the checker exempts FX variables —
//      lo == hi — from the complementarity check entirely).
//   2. Empty column (appears in no active row): fix it at whichever bound
//      minimizes its objective contribution, then remove it.
//   3. Redundant row: a row whose implied activity range — computed from
//      every active column's CURRENT bounds — already lies inside
//      [row_lo, row_hi] regardless of x. Covers the nnz==0 case (implied
//      range is {0}) as well as genuinely-redundant rows with more terms.
//      The row can never bind, so its dual is 0 in every optimal solution —
//      no case analysis needed, y_i = 0 unconditionally.
//   4. Singleton row: a row with exactly one active column j. The row
//      constraint is equivalent to a bound on x_j alone (divide through by
//      the coefficient, flipping the inequality if it's negative); that
//      implied bound is intersected with x_j's own bounds, and the row is
//      removed. Dual recovery (see Postsolve): this row's dual is 0 UNLESS
//      x_j ends up at the specific bound THIS reduction tightened to (and
//      that bound is strictly tighter than what it replaced) — in which
//      case the row is "responsible" for pinning x_j, and its dual is set
//      so x_j's overall reduced cost comes out to exactly 0 (i.e. x_j
//      becomes effectively basic once the row is reinstated).
//   5. Free column singleton: a column with col_lo == -inf and
//      col_hi == +inf that appears in exactly one row, and that row is an
//      EQUALITY (row_lo == row_hi). A free variable is never at a bound in
//      any optimal solution, so it's always basic — its row is really just
//      its own definition in terms of the other columns in that row.
//      Substitute it out (row and column both removed); postsolve
//      recovers its value from the row equation using the other columns'
//      final values, and its row's dual from the zero-reduced-cost
//      condition (y_i = obj_j / a_ij, unconditionally, since a free
//      variable is always basic). Restricted to equality rows specifically
//      so the substituted value is uniquely determined, not merely one of
//      a range of feasible choices.
//   6. General bound tightening (no reduction-stack entry, no postsolve
//      work): for a row with 2+ active columns, each column's implied
//      bound is derived from the row's own bounds and the OTHER active
//      columns' current bounds, and intersected into that column's bounds
//      if strictly tighter. This never removes anything and the row stays
//      a real constraint in the reduced problem, so its dual comes from
//      the ordinary solve — nothing to undo.
//
// Deliberately NOT implemented (documented, not hidden, matching this
// project's existing practice for other partial-scope substitutions):
// forcing rows (would fix multiple variables at once from a single row's
// dual — the general case doesn't reduce to one equation in one unknown
// the way singleton-row and free-column-singleton do), dominated
// column / dual fixing, duplicate row/column detection, and coefficient
// tightening (a MIP-oriented technique with little to say about a
// continuous-only LP, since this project's MILP phase hasn't started).
enum class ReductionType {
  kFixedVariable,
  kEmptyColumn,
  kRedundantRow,
  kSingletonRow,
  kFreeColumnSingleton,
};

struct Reduction {
  ReductionType type;
  int col = -1;   // affected column (all types except kRedundantRow)
  int row = -1;   // affected row (kRedundantRow, kSingletonRow, kFreeColumnSingleton)
  double value = 0.0;  // kFixedVariable/kEmptyColumn: the fixed value.
                        // kFreeColumnSingleton: the row's equality RHS, AS OF
                        // application time — already shifted by every
                        // column fixed earlier in the SAME row (kFixedVariable
                        // folds a fixed column's contribution into row_lo/
                        // row_hi, exactly like it folds into obj_offset).
  double coeff = 0.0;  // a_ij — kSingletonRow, kFreeColumnSingleton.
  // col bounds immediately before/after this reduction — kSingletonRow
  // only, needed by postsolve to tell whether THIS reduction (as opposed
  // to the column's original bound, or an earlier/later one in a chain of
  // stacked tightenings on the same column) is what pinned x_j.
  double prev_lo = 0.0, prev_hi = 0.0;
  double new_lo = 0.0, new_hi = 0.0;
  // kFreeColumnSingleton only: every OTHER column still ACTIVE in this row
  // at application time (excluding `col` itself), i.e. exactly the set
  // `value` was computed against. Postsolve's row equation
  // x_col = (value - sum of these) / coeff must sum over precisely this
  // set — NOT every column the ORIGINAL matrix lists for this row, which
  // would double-count whatever's already folded into `value` (a real bug
  // the full Netlib run caught on greenbeb: a 200+-term row with most of
  // its columns already fixed by the time it became a free-column-singleton
  // candidate).
  std::vector<std::pair<int, double>> other_active_terms;
};

struct PresolveResult {
  core::LpProblem reduced;
  std::vector<Reduction> stack;  // in application order; postsolve undoes in reverse
  bool infeasible = false;       // true if presolve itself proved infeasibility
  std::vector<int> reduced_col_to_original;  // reduced column index -> original column index
  std::vector<int> reduced_row_to_original;  // reduced row index -> original row index
};

PresolveResult Presolve(const core::LpProblem& problem,
                         const core::TolerancePolicy& tol = core::DefaultTolerances());

// Expands a solution computed on `result.reduced` back into `original`'s
// column and row space: fills in every fixed/removed column's known value,
// recomputes row activity directly from the original matrix, and recovers
// a full-length dual vector — reduced-problem rows keep the reduced
// solve's dual value, removed rows get theirs from the rule described
// above for each reduction kind, processed in reverse application order.
core::Solution Postsolve(const core::LpProblem& original, const PresolveResult& result,
                          const core::Solution& reduced_solution,
                          const core::TolerancePolicy& tol = core::DefaultTolerances());

}  // namespace inferno::presolve
