# inferno-solver

SIH 2026 · Problem Statement 26119 · Mangalore Refinery and Petrochemicals
Limited · Team Inferno

A sovereign LP / MILP / QP optimization engine, built from mathematical
foundations — no HiGHS, CBC, GLPK, SCIP, OSQP or SuiteSparse anywhere in the
dependency graph. See [`../BUILD_PLAN_V2.md`](../BUILD_PLAN_V2.md) for the
full 75-day plan, phase gates and kill checkpoints; this README only covers
what is built so far.

## Status: Phase 1.1 — Foundation and proof of life (days 1–3)

What exists:

- `core/` — CSC/CSR sparse matrix, bump-allocator arena, central tolerance policy, the `LpProblem`/`Solution` types every later phase shares
- `io/` — free-form MPS reader (also parses fixed-form Netlib files) and a plain-text solution writer
- `checker/` — the independent solution checker (primal residual, dual residual, complementarity gap), which never shares code with the solver it's checking
- `simplex/dense_simplex.*` — a **throwaway** dense-tableau bounded-variable two-phase simplex, whose only job is to prove the pipeline works end to end. It is O(rows³) per pivot and is expected to be far too slow on the larger Netlib instances; it is replaced by the sparse revised simplex in Phase 2.1 and should not be extended.
- `bench/` — `download_netlib.sh` pulls and decodes the real Netlib LP set from netlib.org; `run_netlib.py` runs the solver over every instance and prints a score

Not yet built (everything else in the phase map): sparse Markowitz LU +
Forrest–Tomlin (Phase 1.2), the real revised simplex + presolve (Phase 2),
PDLP/GPU (Phase 3), MILP + QP (Phase 4), refinery models and packaging
(Phase 5).

## Build

```bash
brew install cmake ninja   # one-time
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Run

```bash
./build/solver bench/netlib/mps/afiro.mps
# AFIRO status=OPTIMAL objective=-464.753 iterations=22 time=0.003s
# checker: PASS primal=... dual=... complementarity=...
```

## Test

```bash
ctest --test-dir build --output-on-failure
```

## Benchmark against the Netlib set

```bash
bash bench/download_netlib.sh   # fetches + decodes into bench/netlib/mps/
python3 bench/run_netlib.py     # runs the solver over all of them, writes bench/results.csv
```

`bench/download_netlib.sh` currently decodes 93 of the ~96 top-level Netlib
instances; `mpc.src`, `stocfor3` and `truss` are shar-bundled Fortran/C
generators on netlib.org rather than plain encoded MPS, and are not yet
wired up (see the comment in that script).

`bench/netlib_optima.csv` only has a verified entry for `afiro` so far —
see that file's header before trusting any pass/fail comparison against
"the published optimum" for other instances; populating it properly is
tracked, not guessed.

## Standing rules (see `BUILD_PLAN_V2.md` for the full list)

1. CI runs the Netlib set on every commit and prints the score.
2. The independent checker (`checker/`) runs on every solve, forever.
3. Every algorithm gets a citation in [`NOTICE_ALGORITHMS.md`](NOTICE_ALGORITHMS.md).
