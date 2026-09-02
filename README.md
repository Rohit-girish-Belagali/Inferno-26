# inferno-solver

SIH 2026 · Problem Statement 26119 · Mangalore Refinery and Petrochemicals
Limited · Team Inferno

A sovereign LP / MILP / QP optimization engine, built from mathematical
foundations — no HiGHS, CBC, GLPK, SCIP, OSQP or SuiteSparse anywhere in the
dependency graph. See [`../BUILD_PLAN_V2.md`](../BUILD_PLAN_V2.md) for the
full 75-day plan, phase gates and kill checkpoints; this README only covers
what is built so far.

## Status: Phase 2 — revised simplex + presolve (days 11–32)

What exists:

- `core/` — CSC/CSR sparse matrix, bump-allocator arena, central tolerance policy, the `LpProblem`/`Solution` types every later phase shares
- `io/` — MPS reader: tries free-form (whitespace-tokenized) first, falls back to strict fixed-column parsing for the minority of older Netlib files that need it (embedded spaces in names, blank continuation fields) — plus a plain-text solution writer
- `checker/` — the independent solution checker (primal residual, dual residual, complementarity gap — fixed variables correctly exempted from the complementarity condition), which never shares code with the solver it's checking
- `la/` — the linear algebra spine: geometric-mean scaling, sparse Markowitz LU with threshold pivoting, Gilbert-Peierls FTRAN/BTRAN, and a basis-update path. **The update is product-form-of-the-inverse (PFI), not full Forrest-Tomlin** — see `la/basis_factorization.hpp`'s header comment for why. Paired with a refactorization policy that bounds the eta chain.
- `simplex/revised_simplex.*` — **the real solver now**: bounded-variable primal simplex built on `la/`'s sparse LU instead of a dense tableau, with scaling applied around the whole solve. A greedy weighted-matching crash basis (not the all-slack start), Devex pricing (with a safe deterministic tie-breaking perturbation, and a Bland's-rule fallback after sustained degenerate pivots), a two-pass ratio test (minimum step, then most numerically stable among ties — not yet full Harris). Solves **80/93 Netlib instances with `--presolve`, zero checker-rejected "optimal" claims** — see `NOTICE_ALGORITHMS.md` for the honest, five-increment measurement (Dantzig → Devex → +scaling → +crash basis → +perturbed tie-breaking) that got there.
- `simplex/dense_simplex.*` — the Phase 1.1 **throwaway** dense tableau, kept only as a comparison baseline (`--solver dense`); no longer the default and not extended further.
- `presolve/presolve.*` — six reduction kinds to a fixpoint with exact postsolve: fixed-variable substitution, empty-column removal, redundant-row removal, singleton-row bound-tightening-then-removal, free-column-singleton substitution, and general (non-removing) bound tightening — the last three needed real postsolve dual (and, for free column singletons, primal) recovery, not just bound narrowing, see `NOTICE_ALGORITHMS.md` for the derivations and the real bugs the full Netlib run caught fixing them. Verified via the plan's own standing rule (presolve on/off equivalence, checked against real Netlib data) — `--presolve` is opt-in, not yet the CLI/bench default, since the reduction set is still incomplete (forcing rows, dominated columns, duplicate rows/columns, coefficient tightening remain; duplicate-column merging was attempted and reverted after a deep cross-reduction ordering bug, also in `NOTICE_ALGORITHMS.md`).
- `simplex/dual_simplex.*` — bounded-variable dual simplex on the same LU spine, for BUILD_PLAN_V2.md's "not optional" (MILP node warm starts, Phase 4). **Only problems with a trivial dual-feasible start are supported** — a general dual phase 1 isn't implemented, see `NOTICE_ALGORITHMS.md`. 35/93 verified, zero checker-rejected claims (`--solver dual`).
- `dashboard/` — two linked pages, `index.html` (wireframe module 6, a live-demo console walking the real Prepare → Solve → Verify pipeline; the GPU/PDLP lane is honestly labeled simulated — Phase 3 hasn't been built, no GPU hardware on the dev machine) and `report.html` (module 5, a sortable/filterable view of the complete real `bench/results.csv`, all 93 instances). Both use genuine numbers from this project's own bench runs, never invented.
- `bench/` — `download_netlib.sh` pulls and decodes the real Netlib LP set from netlib.org; `run_netlib.py --solver {revised,dense,dual,managed} [--presolve]` runs a solve path over every instance and prints a score, scoring by independent-checker PASS rather than solver-claimed status

Not yet built: steepest-edge pricing, full Harris ratio test, a
bound-flipping dual ratio test, a general dual phase 1, forcing
rows/dominated columns/duplicate detection/coefficient tightening for
presolve (rest of Phase 2); PDLP/GPU (Phase 3); MILP + QP (Phase 4);
refinery models and packaging (Phase 5). True Forrest-Tomlin (upgrading
from the current PFI update) is also outstanding, tracked in
`NOTICE_ALGORITHMS.md`.

## Build

```bash
brew install cmake ninja   # one-time
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Run

```bash
./build/solver bench/netlib/mps/afiro.mps
# AFIRO solver=revised status=OPTIMAL objective=-464.753 iterations=27 time=0.0004s
# checker: PASS primal=... dual=... complementarity=...

./build/solver bench/netlib/mps/afiro.mps --solver dense   # Phase 1.1 baseline, for comparison
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

## Linear algebra (`la/`)

`ctest` covers `la/` with hand-built, identity, singular, threshold-pivoting
and random 60×60 sparse matrices, plus a 1000-update stress test checked
against the Phase 1 gate (updates match a fresh refactorization within
1e-8, and run ≥20x faster — currently ~44.5x). It has also been checked
against **real** optimal bases: solve a Netlib instance with the Phase 1.1
dense simplex, reconstruct the actual basis matrix it landed on, and
confirm the new sparse LU factors and solves it too (residuals at machine
epsilon). That real-basis check isn't part of `ctest` (it needs the
downloaded Netlib set) — see `la/manual_checks/`.

## Standing rules (see `BUILD_PLAN_V2.md` for the full list)

1. CI runs the Netlib set on every commit and prints the score.
2. The independent checker (`checker/`) runs on every solve, forever.
3. Every algorithm gets a citation in [`NOTICE_ALGORITHMS.md`](NOTICE_ALGORITHMS.md).
