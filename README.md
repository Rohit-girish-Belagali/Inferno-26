# inferno-solver

SIH 2026 · Problem Statement 26119 · Mangalore Refinery and Petrochemicals
Limited · Team Inferno

A sovereign LP / MILP / QP optimization engine, built from mathematical
foundations — no HiGHS, CBC, GLPK, SCIP, OSQP or SuiteSparse anywhere in the
dependency graph. See [`../BUILD_PLAN_V2.md`](../BUILD_PLAN_V2.md) for the
full 75-day plan, phase gates and kill checkpoints; this README only covers
what is built so far.

## Status: Phase 2.1 — the real revised simplex (days 11–25)

What exists:

- `core/` — CSC/CSR sparse matrix, bump-allocator arena, central tolerance policy, the `LpProblem`/`Solution` types every later phase shares
- `io/` — MPS reader: tries free-form (whitespace-tokenized) first, falls back to strict fixed-column parsing for the minority of older Netlib files that need it (embedded spaces in names, blank continuation fields) — plus a plain-text solution writer
- `checker/` — the independent solution checker (primal residual, dual residual, complementarity gap — fixed variables correctly exempted from the complementarity condition), which never shares code with the solver it's checking
- `la/` — the linear algebra spine: geometric-mean scaling, sparse Markowitz LU with threshold pivoting, Gilbert-Peierls FTRAN/BTRAN, and a basis-update path. **The update is product-form-of-the-inverse (PFI), not full Forrest-Tomlin** — see `la/basis_factorization.hpp`'s header comment for why. Paired with a refactorization policy that bounds the eta chain.
- `simplex/revised_simplex.*` — **the real solver now**: bounded-variable primal simplex built on `la/`'s sparse LU instead of a dense tableau. Dantzig pricing (Bland's-rule fallback after sustained degenerate pivots), a two-pass ratio test (minimum step, then most numerically stable among ties — not yet full Harris). See `NOTICE_ALGORITHMS.md`.
- `simplex/dense_simplex.*` — the Phase 1.1 **throwaway** dense tableau, kept only as a comparison baseline (`--solver dense`); no longer the default and not extended further.
- `dashboard/` — an interactive demo instrument (`dashboard/index.html`) walking the real Prepare → Solve → Verify pipeline against solved instances, using genuine numbers from this project's own bench runs; the GPU/PDLP lane is honestly labeled simulated (Phase 3 hasn't been built — no GPU hardware on the dev machine).
- `bench/` — `download_netlib.sh` pulls and decodes the real Netlib LP set from netlib.org; `run_netlib.py --solver {revised,dense}` runs either solve path over every instance and prints a score, scoring by independent-checker PASS rather than solver-claimed status

Not yet built: presolve, Devex/steepest-edge pricing, full Harris ratio
test, dual simplex (rest of Phase 2); PDLP/GPU (Phase 3); MILP + QP
(Phase 4); refinery models and packaging (Phase 5). True Forrest-Tomlin
(upgrading from the current PFI update) is also outstanding, tracked in
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
# AFIRO solver=revised status=OPTIMAL objective=-464.753 iterations=18 time=0.0002s
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
