/* Inferno solver — stable C ABI (BUILD_PLAN_V2.md Phase 5).
 *
 * A C interface rather than the C++ one for the usual reason: C++ has no
 * stable ABI across compilers or standard-library versions, so anything
 * that wants to link this engine from another language, or from a binary
 * built by a different toolchain, needs a C surface. This is that surface.
 *
 * Conventions, all deliberate:
 *   - Every function that can fail returns an inferno_status. Nothing
 *     throws across this boundary; C++ exceptions are caught at the edge
 *     and converted, because letting one unwind into a C caller is
 *     undefined behaviour.
 *   - Ownership is explicit: whatever you create, you destroy. No function
 *     here returns memory the caller must free -- solution data is copied
 *     into caller-provided buffers instead, so there is no allocator
 *     mismatch across a library boundary.
 *   - Indices are 0-based. Bounds use INFERNO_INFINITY for "unbounded".
 *   - All pointers must be non-NULL unless documented otherwise; passing
 *     NULL returns INFERNO_STATUS_INVALID_ARGUMENT rather than crashing.
 */
#ifndef INFERNO_H
#define INFERNO_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define INFERNO_INFINITY (1.0e308)

typedef enum {
  INFERNO_STATUS_OK = 0,
  INFERNO_STATUS_INVALID_ARGUMENT = 1,
  INFERNO_STATUS_READ_FAILED = 2,
  INFERNO_STATUS_NOT_SOLVED = 3,
  INFERNO_STATUS_INTERNAL_ERROR = 4
} inferno_status;

/* The outcome of a solve. Distinct from inferno_status: a solve can
 * succeed as an API call while proving the problem infeasible. */
typedef enum {
  INFERNO_SOLVE_OPTIMAL = 0,
  INFERNO_SOLVE_INFEASIBLE = 1,
  INFERNO_SOLVE_UNBOUNDED = 2,
  INFERNO_SOLVE_ITERATION_LIMIT = 3,
  INFERNO_SOLVE_NUMERICAL_ERROR = 4
} inferno_solve_status;

typedef struct inferno_problem inferno_problem;

/* --- Lifecycle --- */
inferno_problem* inferno_problem_create(void);
void inferno_problem_destroy(inferno_problem* p);

/* Loads an MPS file. Returns INFERNO_STATUS_READ_FAILED if the file cannot
 * be parsed; the problem is left unmodified in that case. */
inferno_status inferno_problem_load_mps(inferno_problem* p, const char* path);

/* Builds a problem directly, for callers generating models in memory.
 * The matrix is compressed sparse column: col_ptr has num_cols + 1 entries,
 * row_idx and values have col_ptr[num_cols] entries. Everything is copied,
 * so the caller's buffers may be freed immediately after this returns. */
inferno_status inferno_problem_set(inferno_problem* p, int num_rows, int num_cols,
                                    const int* col_ptr, const int* row_idx,
                                    const double* values, const double* obj,
                                    const double* col_lo, const double* col_hi,
                                    const double* row_lo, const double* row_hi);

/* --- Solving --- */
/* Solves with the revised simplex. `presolve` non-zero enables the
 * presolve/postsolve pass. The independent checker always runs; its
 * verdict is available via inferno_get_checker_passed, and a caller that
 * cares about correctness should consult it rather than trusting the
 * solve status alone -- that is the whole discipline of this project. */
inferno_status inferno_solve(inferno_problem* p, int presolve);

inferno_status inferno_get_solve_status(const inferno_problem* p, inferno_solve_status* out);
inferno_status inferno_get_objective(const inferno_problem* p, double* out);
inferno_status inferno_get_iterations(const inferno_problem* p, int* out);
/* Non-zero if the independent checker accepted the solution. */
inferno_status inferno_get_checker_passed(const inferno_problem* p, int* out);

/* Copies the primal solution into `out`, which must have room for
 * `count` doubles; `count` must equal num_cols. */
inferno_status inferno_get_solution(const inferno_problem* p, double* out, int count);
/* Copies the row duals (shadow prices); `count` must equal num_rows. */
inferno_status inferno_get_duals(const inferno_problem* p, double* out, int count);

inferno_status inferno_get_num_rows(const inferno_problem* p, int* out);
inferno_status inferno_get_num_cols(const inferno_problem* p, int* out);

/* Static, human-readable description of a status code. Never NULL, and
 * points to a string literal the caller must not free. */
const char* inferno_status_string(inferno_status s);

/* Version, for callers that need to gate on capabilities. */
const char* inferno_version(void);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* INFERNO_H */
