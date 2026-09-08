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


  /* --- MILP through the C ABI, plus the node-event interface. ---------
   * Knapsack: maximize 5a + 4b (as minimize -5a - 4b) subject to
   * 6a + 4b <= 24 with a, b integer in [0, 3]. The optimum is a=2, b=3
   * giving 22, so the minimized objective is -22. */
  {
    inferno_problem* mp = inferno_problem_create();
    expect(mp != NULL, "create mip problem");

    int col_ptr[3] = {0, 1, 2};
    int row_idx[2] = {0, 0};
    double values[2] = {6.0, 4.0};
    double obj[2] = {-5.0, -4.0};
    double col_lo[2] = {0.0, 0.0};
    double col_hi[2] = {3.0, 3.0};
    double row_lo[1] = {-INFERNO_INFINITY};
    double row_hi[1] = {24.0};
    int flags[2] = {1, 1};

    expect(inferno_problem_set(mp, 1, 2, col_ptr, row_idx, values, obj, col_lo, col_hi,
                              row_lo, row_hi) == INFERNO_STATUS_OK,
          "set mip problem");

    /* MILP accessors must refuse to answer before a MILP solve has run. */
    double tmp = 0.0;
    expect(inferno_get_best_bound(mp, &tmp) == INFERNO_STATUS_NOT_SOLVED,
          "best bound is unavailable before solving");

    expect(inferno_problem_set_integer(mp, flags, 2) == INFERNO_STATUS_OK, "set integer flags");
    expect(inferno_problem_set_integer(mp, flags, 5) == INFERNO_STATUS_INVALID_ARGUMENT,
          "wrong integer flag count is rejected");

    expect(inferno_solve_mip(mp, 10.0, 100000, 0.0, 1) == INFERNO_STATUS_OK, "solve mip");

    inferno_solve_status st;
    expect(inferno_get_solve_status(mp, &st) == INFERNO_STATUS_OK, "mip status");
    expect(st == INFERNO_SOLVE_OPTIMAL, "mip is optimal");

    double objective = 0.0, bound = 0.0, gap = 1.0;
    int nodes = -1, proved = 0, passed = 0;
    expect(inferno_get_objective(mp, &objective) == INFERNO_STATUS_OK, "mip objective");
    expect(fabs(objective - (-22.0)) < 1e-9, "mip objective is -22");
    expect(inferno_get_best_bound(mp, &bound) == INFERNO_STATUS_OK, "mip bound");
    expect(bound <= objective + 1e-9, "bound does not exceed the objective");
    expect(inferno_get_gap(mp, &gap) == INFERNO_STATUS_OK, "mip gap");
    expect(gap < 1e-6, "gap is closed");
    expect(inferno_get_nodes_explored(mp, &nodes) == INFERNO_STATUS_OK, "node count");
    expect(nodes >= 0, "node count is non-negative");
    expect(inferno_get_proved_optimal(mp, &proved) == INFERNO_STATUS_OK, "proved flag");
    expect(proved == 1, "optimality was proved");
    expect(inferno_get_checker_passed(mp, &passed) == INFERNO_STATUS_OK, "mip checker");
    expect(passed == 1, "independent MILP checker passed");

    double xs[2] = {0.0, 0.0};
    expect(inferno_get_solution(mp, xs, 2) == INFERNO_STATUS_OK, "mip solution");
    expect(fabs(xs[0] - 2.0) < 1e-9 && fabs(xs[1] - 3.0) < 1e-9, "mip solution is (2, 3)");

    /* The event interface: at least the root must have been reported, and
     * every event must carry a known outcome. */
    int count = -1;
    expect(inferno_get_event_count(mp, &count) == INFERNO_STATUS_OK, "event count");
    expect(count >= 1, "at least the root relaxation was reported");

    inferno_node_event evs[64];
    int written = -1;
    expect(inferno_get_events(mp, evs, 0, 64, &written) == INFERNO_STATUS_OK, "read events");
    expect(written >= 1 && written <= count, "event read count is sane");
    expect(evs[0].outcome == INFERNO_NODE_ROOT, "the first event is the root relaxation");
    expect(evs[0].depth == 0, "the root is at depth 0");
    for (int i = 0; i < written; ++i) {
      expect(evs[i].outcome >= INFERNO_NODE_ROOT && evs[i].outcome <= INFERNO_NODE_RELAXATION_FAILED,
            "event outcome is a known code");
      expect(evs[i].depth >= 0, "event depth is non-negative");
      expect(inferno_node_outcome_string(evs[i].outcome) != NULL, "outcome has a name");
    }

    /* Reading past the end must report zero written, not walk off the array. */
    expect(inferno_get_events(mp, evs, count + 10, 64, &written) == INFERNO_STATUS_OK,
          "read beyond the end is not an error");
    expect(written == 0, "reading beyond the end writes nothing");

    inferno_problem_destroy(mp);
  }

  if (g_failures) { fprintf(stderr, "%d test(s) failed\n", g_failures); return 1; }
  printf("all C ABI tests passed\n");
  return 0;
}
