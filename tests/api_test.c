/* C ABI test — deliberately compiled as C, not C++, because the entire
 * point of the ABI is that a C consumer can use it. If this ever stops
 * compiling as C, the header has grown a C++-ism and the ABI is broken. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "api/inferno.h"

static int g_failures = 0;

static void expect(int cond, const char* msg) {
  if (!cond) { fprintf(stderr, "FAIL: %s\n", msg); ++g_failures; }
  else { printf("ok: %s\n", msg); }
}

int main(void) {
  printf("inferno version %s\n", inferno_version());

  /* --- Failure paths first: an API is judged by what it does when
   * misused, and none of these may crash. --- */
  expect(inferno_problem_load_mps(NULL, "x") == INFERNO_STATUS_INVALID_ARGUMENT,
         "NULL problem handle is rejected, not dereferenced");
  {
    inferno_problem* p = inferno_problem_create();
    expect(p != NULL, "problem creates");
    expect(inferno_problem_load_mps(p, "/definitely/not/a/file.mps") == INFERNO_STATUS_READ_FAILED,
           "missing file reports READ_FAILED");
    double obj;
    expect(inferno_get_objective(p, &obj) == INFERNO_STATUS_NOT_SOLVED,
           "querying before solving reports NOT_SOLVED");
    inferno_problem_destroy(p);
  }
  inferno_problem_destroy(NULL);  /* must be a no-op, like free(NULL) */
  printf("ok: destroying NULL is a no-op\n");

  /* --- A real solve, built through the in-memory path.
   *   max 3x + 2y   ->   min -3x - 2y
   *   s.t.  x +  y <= 4
   *         x + 3y <= 6
   *         0 <= x,y <= 3
   * The optimum is x=3, y=1, objective -11. --- */
  {
    inferno_problem* p = inferno_problem_create();
    int    col_ptr[3] = {0, 2, 4};
    int    row_idx[4] = {0, 1, 0, 1};
    double values[4]  = {1.0, 1.0, 1.0, 3.0};
    double obj[2]     = {-3.0, -2.0};
    double col_lo[2]  = {0.0, 0.0};
    double col_hi[2]  = {3.0, 3.0};
    double row_lo[2]  = {-INFERNO_INFINITY, -INFERNO_INFINITY};
    double row_hi[2]  = {4.0, 6.0};

    inferno_status st = inferno_problem_set(p, 2, 2, col_ptr, row_idx, values, obj,
                                             col_lo, col_hi, row_lo, row_hi);
    expect(st == INFERNO_STATUS_OK, "problem_set accepts a well-formed model");

    int nr = 0, nc = 0;
    inferno_get_num_rows(p, &nr);
    inferno_get_num_cols(p, &nc);
    expect(nr == 2 && nc == 2, "dimensions round-trip");

    expect(inferno_solve(p, 0) == INFERNO_STATUS_OK, "solve returns OK");

    inferno_solve_status ss;
    expect(inferno_get_solve_status(p, &ss) == INFERNO_STATUS_OK && ss == INFERNO_SOLVE_OPTIMAL,
           "solve reports OPTIMAL");

    int passed = 0;
    inferno_get_checker_passed(p, &passed);
    expect(passed == 1, "independent checker accepted the solution through the ABI");

    double v = 0.0;
    inferno_get_objective(p, &v);
    expect(fabs(v - (-11.0)) < 1e-6, "objective is -11");

    double x[2] = {0.0, 0.0};
    expect(inferno_get_solution(p, x, 2) == INFERNO_STATUS_OK, "solution copies out");
    expect(fabs(x[0] - 3.0) < 1e-6 && fabs(x[1] - 1.0) < 1e-6, "solution is x=3, y=1");

    expect(inferno_get_solution(p, x, 5) == INFERNO_STATUS_INVALID_ARGUMENT,
           "a wrong buffer size is rejected rather than overrunning");

    double y[2];
    expect(inferno_get_duals(p, y, 2) == INFERNO_STATUS_OK, "duals copy out");

    inferno_problem_destroy(p);
  }

  /* --- A malformed matrix must be rejected without corrupting state. --- */
  {
    inferno_problem* p = inferno_problem_create();
    int    col_ptr[2] = {0, 1};
    int    row_idx[1] = {99};       /* row 99 does not exist */
    double values[1]  = {1.0};
    double obj[1]     = {1.0};
    double lo[1] = {0.0}, hi[1] = {1.0};
    double rlo[1] = {0.0}, rhi[1] = {1.0};
    expect(inferno_problem_set(p, 1, 1, col_ptr, row_idx, values, obj, lo, hi, rlo, rhi)
               == INFERNO_STATUS_INVALID_ARGUMENT,
           "out-of-range row index is rejected");
    inferno_problem_destroy(p);
  }

  expect(strcmp(inferno_status_string(INFERNO_STATUS_OK), "ok") == 0, "status strings work");

  if (g_failures) { fprintf(stderr, "%d test(s) failed\n", g_failures); return 1; }
  printf("all C ABI tests passed\n");
  return 0;
}
