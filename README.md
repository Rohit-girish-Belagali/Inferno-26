# inferno-solver

SIH 2026 · Problem Statement 26119 · Mangalore Refinery and Petrochemicals
Limited · Team Inferno

A sovereign LP / MILP / QP optimization engine, built from mathematical
foundations — no HiGHS, CBC, GLPK, SCIP, OSQP or SuiteSparse anywhere in the
dependency graph. See [`../BUILD_PLAN_V2.md`](../BUILD_PLAN_V2.md) for the
full 75-day plan, phase gates and kill checkpoints; this README only covers
what is built so far. [`ARCHITECTURE.md`](ARCHITECTURE.md) explains how the
pieces fit together and why each is shaped the way it is.

## Status: LP 91/93 Netlib · QP via ADMM · PDLP CPU reference · refinery model · C ABI + Python

What exists:

- `core/` — CSC/CSR sparse matrix, bump-allocator arena, central tolerance policy, the `LpProblem`/`Solution` types every later phase shares
- `io/` — MPS reader: tries free-form (whitespace-tokenized) first, falls back to strict fixed-column parsing for the minority of older Netlib files that need it (embedded spaces in names, blank continuation fields) — plus a plain-text solution writer
- `checker/` — the independent solution checker (primal residual, dual residual, complementarity gap — fixed variables correctly exempted from the complementarity condition), which never shares code with the solver it's checking
- `la/` — the linear algebra spine: geometric-mean scaling, sparse Markowitz LU with threshold pivoting, Gilbert-Peierls FTRAN/BTRAN, and a basis-update path. **The update is product-form-of-the-inverse (PFI), not full Forrest-Tomlin** — see `la/basis_factorization.hpp`'s header comment for why. Paired with a refactorization policy that bounds the eta chain.
- `simplex/revised_simplex.*` — **the real solver now**: bounded-variable primal simplex built on `la/`'s sparse LU instead of a dense tableau, with scaling applied around the whole solve. A greedy weighted-matching crash basis (not the all-slack start), Devex pricing (with a safe deterministic tie-breaking perturbation, and a Bland's-rule fallback after sustained degenerate pivots), a Harris two-pass ratio test with adaptive tolerance expansion (which is what clears degenerate stalling), and singular-basis repair: when a mid-solve refactorization finds the basis singular, the LU reports which columns it could not pivot and those are swapped for the slacks of the uncovered rows (provably nonsingular) so the solve continues instead of being abandoned. Solves **91/93 Netlib instances with `--presolve`, zero checker-rejected "optimal" claims, zero numerical errors, and no wrong answers** — the two remaining (`dfl001`, `greenbea`) are wall-clock timeouts on the largest models, not correctness failures. See `NOTICE_ALGORITHMS.md` for the honest, increment-by-increment measurement that got there (Dantzig → Devex → +scaling → +crash basis → +perturbed tie-breaking → +singular-basis repair → +LU pivot search and refactorization interval → +Harris expansion), including the two changes that were measured, found to cost an instance, and reverted rather than shipped. All 91 verified instances solve in 61.3s total. What separates this from the plan's 98/98 gate is now raw throughput on the two largest models — steepest-edge pricing is the next lever.
- `simplex/dense_simplex.*` — the Phase 1.1 **throwaway** dense tableau, kept only as a comparison baseline (`--solver dense`); no longer the default and not extended further.
- `presolve/presolve.*` — six reduction kinds to a fixpoint with exact postsolve: fixed-variable substitution, empty-column removal, redundant-row removal, singleton-row bound-tightening-then-removal, free-column-singleton substitution, and general (non-removing) bound tightening — the last three needed real postsolve dual (and, for free column singletons, primal) recovery, not just bound narrowing, see `NOTICE_ALGORITHMS.md` for the derivations and the real bugs the full Netlib run caught fixing them. Verified via the plan's own standing rule (presolve on/off equivalence, checked against real Netlib data) — `--presolve` is opt-in, not yet the CLI/bench default, since the reduction set is still incomplete (forcing rows, dominated columns, duplicate rows/columns, coefficient tightening remain; duplicate-column merging was attempted and reverted after a deep cross-reduction ordering bug, also in `NOTICE_ALGORITHMS.md`).
- `simplex/dual_simplex.*` — bounded-variable dual simplex on the same LU spine, for BUILD_PLAN_V2.md's "not optional" (MILP node warm starts, Phase 4). **Only problems with a trivial dual-feasible start are supported** — a general dual phase 1 isn't implemented, see `NOTICE_ALGORITHMS.md`. 35/93 verified, zero checker-rejected claims (`--solver dual`).
- `dashboard/` — two linked pages, `index.html` (wireframe module 6, a live-demo console walking the real Prepare → Solve → Verify pipeline; the GPU lane is honestly labeled simulated — the PDLP CPU reference is real, but there is no GPU on this machine and no GPU number is claimed) and `report.html` (module 5, a sortable/filterable view of the complete real `bench/results.csv`, all 93 instances). Both use genuine numbers from this project's own bench runs, never invented.
- `firstorder/pdlp.*` — **PDLP**, a primal-dual first-order solver (Phase 3's CPU reference). PDHG on the saddle-point form, the dual prox reduced to a clamp by Moreau decomposition, power-iteration norm estimate, averaged iterates, adaptive restarts, primal weight balancing, and preconditioning. Measured at **4 of 15** small Netlib instances checker-verified within 30k iterations (matching the simplex objective to ~1e-10); the rest stall between 1e-5 and 1e-1 relative KKT, which needs the adaptive step size and proper restart criterion that are not built yet. A reference implementation, not a competitive LP solver. **The CUDA port is hardware-blocked** (no NVIDIA GPU here) and no GPU speedup is claimed — see `STATUS.md`.
- `models/breadth_models.*` — **transportation and economic dispatch** (Phase 5 breadth models), showing the engine is a general LP engine rather than something shaped around one use case. The dispatch model is named honestly: unit commitment proper is a MILP, so this is the continuous economic-dispatch sub-problem, with ramp limits coupling consecutive periods.
- `models/crude_blending.*` + `src/refinery_demo.cpp` — the **refinery crude-blending model** (Phase 5). An LP: the volume-weighted quality ratio is cleared to a linear constraint, which is what keeps it solvable without MILP. `./build/refinery_demo` prints the blend table, realised properties against spec, and marginal values in profit terms.
- `qp/admm.*` — **convex QP via ADMM** (Phase 4Q), completing the LP/QP pair the problem statement names. Splits `z = Ax` so the quadratic and the inequalities are never handled together; the x-step KKT matrix is quasi-definite and fixed across iterations, so it is factorized once with the Phase 1.2 LU and every iteration is two triangular solves plus vector work. Tested against hand-derived optima.
- `bench/` — `download_netlib.sh` pulls and decodes the real Netlib LP set from netlib.org; `run_netlib.py --solver {revised,dense,dual,managed} [--presolve]` runs a solve path over every instance and prints a score, scoring by independent-checker PASS rather than solver-claimed status

Not yet built, stated plainly: steepest-edge pricing (attempted and
reverted — see `NOTICE_ALGORITHMS.md`), a bound-flipping dual ratio test,
a general dual phase 1, and the remaining presolve reductions (forcing
rows, dominated columns, duplicate detection, coefficient tightening).
PDLP's adaptive step size and restart criterion, and the CUDA port, which
is hardware-blocked. **MILP (Phase 4) is deliberately out of scope** — the plan's own day-25 kill checkpoint says to abandon MILP
unless all 98 Netlib instances solve, and at 91/93 they do not; see
`STATUS.md`. From Phase 5: the unit-scheduling MILP and a
Dolan–Moré profile *against HiGHS* (the one here is a runtime
distribution, since with a single solver the reference time is its own
best). True
Forrest–Tomlin (upgrading from the current PFI update) is also
outstanding, tracked in `NOTICE_ALGORITHMS.md`.

## The demo

```bash
./demo/run_demo.sh
```

Ten steps, walking the analyst user flow exactly as the PS 26119 proposal
deck describes it — model in (MPS) → parse & validate → presolve → scale →
concurrent solve manager → exact postsolve → independent checker → report
— then the refinery model, the breadth models, QP, all three modeling APIs,
the full benchmark, and the clean-room proof.

It ends with a **capability scorecard against the PS's six named modules,
including the ones that are not built**. The GPU acceleration layer is the
headline item in the problem statement's own title and it is marked NOT
BUILT, because there is no NVIDIA GPU on this machine and no speedup is
claimed. A demo that hides its gaps is worth less than one that states
them — the first question from any panel is the gap.

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

./build/refinery_demo    # the Phase 5 crude-blending model, solved and verified
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

### Python

```python
import inferno                     # python/ on your PYTHONPATH
p = inferno.Problem.from_mps("bench/netlib/mps/afiro.mps")
r = p.solve(presolve=True)
print(r.status_name, r.objective, r.checker_passed)
# OPTIMAL -464.75314285714353 True
```

Bindings are `ctypes` over the C ABI, not a compiled extension — no build
step, no compiler needed on the user's machine, and no matching of
Python's ABI, so it works against a prebuilt `libinferno` on any Python 3.
`r.checker_passed` is the field to trust: it carries the independent
checker's verdict, not the solver's own claim.

### Reproducing the headline numbers

Everything this README claims is regenerated by these four commands — the
report reads only the committed `bench/results.csv` and computes each
figure rather than transcribing it, so if a number here and one it prints
ever disagree, the report is right.

```bash
bash bench/download_netlib.sh                              # fetch + decode the real Netlib set
python3 bench/run_netlib.py --solver revised --presolve    # solve all 93, write results.csv
python3 bench/make_report.py                               # the benchmark report
bash bench/verify_clean_room.sh                            # prove no third-party solver is linked
```

`verify_clean_room.sh` is the mechanical form of this project's central
claim. It greps every source and build file for the name of any
third-party solver, inspects what the built binaries actually link with
`otool`/`ldd`, and lists every CMake link dependency. It exits non-zero if
anything is found, so the claim is checkable rather than asserted.

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
