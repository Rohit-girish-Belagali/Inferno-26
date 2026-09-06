# Architecture

How the solver is put together, and — more usefully — why each piece is
the way it is. Every design decision here was made against a measurement
or a correctness argument, and where a decision was reversed, the reversal
is recorded too. `NOTICE_ALGORITHMS.md` carries the citations and the
per-component measurement history; this document is the shape of the
thing.

## The one idea the whole project is built around

**The component that produces an answer never gets to decide whether the
answer is correct.**

`checker/` recomputes the primal residual, the dual residual and the
complementarity gap from the raw problem data. It shares no code with any
solver. It runs on every solve, in every mode, including through the C
ABI — and the benchmark harness scores by the checker's verdict, not by
what the solver claimed. An instance the solver calls `OPTIMAL` and the
checker rejects is counted as a **wrong answer**, not folded into the
success column.

This is not defensive decoration. Over this project's history it caught:

- three instances where the dual simplex reported a confident optimum with
  large complementarity violations, because dual feasibility was assumed
  to be an invariant and never re-verified;
- a presolve postsolve bug that made `greenbeb`'s recovered solution
  primal-infeasible by ~10⁴;
- an objective-substitution term dropped for free column singletons.

None of these were visible from the solver's own point of view. All of
them looked like success.

## Layers

```
        io/            models/          api/  python/
     MPS reader     refinery, breadth   C ABI, ctypes
          \              |               /
           \             |              /
            +---------- core/ ---------+        checker/
              LpProblem, QpProblem,              independent
              CSC/CSR sparse, tolerances         verification
                         |
        +----------------+-----------------+
        |                |                 |
    presolve/        simplex/          firstorder/   qp/
    6 reductions   revised, dual,        PDLP        ADMM
    + postsolve    solve manager
                         |
                        la/
        scaling · Markowitz LU · Gilbert-Peierls
        FTRAN/BTRAN · PFI basis update
```

`la/` is the spine. Every higher layer that needs to solve a linear system
goes through it — simplex for FTRAN/BTRAN, and QP's ADMM for its KKT
factorization. That reuse was a plan decision ("the KKT factorization from
Phase 1.2 is reusable"), and it held: the QP solver is thin precisely
because the hard part already existed.

## Why the pieces are shaped this way

**Sparse from the bottom up.** `core/` stores compressed sparse column
because FTRAN, BTRAN and pivoting are all column operations. A row-major
view is built alongside it where a genuinely row-oriented question is
asked — the Devex weight update needs the pivot *row* against every
nonbasic column, and answering that column-by-column against CSC meant
touching every nonzero in the problem on every pivot.

**PFI instead of Forrest–Tomlin.** The plan names Forrest–Tomlin's LU-form
update. Its exact permutation-and-elimination procedure was judged too
easy to get subtly wrong from memory without a reference to check against,
so the simpler product-form update was implemented instead and paired with
a refactorization policy that bounds the eta chain PFI alone does not
guarantee against. This is a documented substitution, not an oversight,
and it is still the honest state of the code.

**The refactorization interval scales with basis size.** Component timing
found refactorization was **77%** of `pilot87`'s entire solve. A flat cap
is the wrong shape: a rebuild costs more as *m* grows, so a bigger basis
should amortize it over more updates. Measured at 100/200/400/800, the
curve is sharp in both directions — 800 did not finish inside ten minutes,
because past some length the eta chain makes every solve against the basis
ruinous. "Just refactorize less" stops being true well before chain growth
does.

**The LU pivot search and that interval are one decision, not two.**
Bucketing columns by nonzero count made factorization ~5× cheaper but
*cost* an instance on its own. The reason was not the search. A cheaper
rebuild changes what the right rebuild *interval* is — with a cheap
factorization you rebuild more often and keep the chain short, which
speeds up every solve too. Landed together they gained what neither gained
alone.

**Degeneracy is handled by an expanding tolerance, and that expansion is
strictly non-losing.** The Harris ratio test measures against bounds
relaxed by a delta, so a variable already sitting on its bound still
admits a positive step — the whole degenerate case. A flat delta is wrong
at both ends, measured: too small and `tuff` burns its entire iteration
budget; too large and `degen3` fails numerically. So it ramps only while
degenerate pivots accumulate. And because at one ramp setting `degen3`
reported `UNBOUNDED` — a wrong answer the checker *cannot* catch, since it
only validates claimed optima — any attempt ending in error, unboundedness
or an iteration limit redoes the whole solve from scratch with the
expansion off. From scratch matters: an in-place restart carries the
damaged basis forward.

**A singular basis is repaired, not surrendered to.** When refactorization
finds the basis singular, the LU reports which columns it could not pivot
and which rows were left uncovered; swapping those columns for the slacks
of the uncovered rows is provably nonsingular. Abandoning thousands of
good iterations because one basis went singular was never the right
response.

## What is deliberately absent

- **MILP.** The plan's own day-25 kill checkpoint: *"do all 98 Netlib
  instances solve to 1e-6? No → abandon MILP."* At 91/93 it does not, so
  MILP is unbuilt. Its stated rationale — a correct LP solver with an
  honest scope beats a broken MILP solver with an ambitious one — is
  followed rather than overridden.
- **A GPU number.** The CUDA port is hardware-blocked (no NVIDIA GPU on
  the development machine). The day-40 checkpoint says drop the claim
  rather than fabricate a speedup, so there is no GPU figure anywhere,
  real or claimed.
- **A comparison against HiGHS.** The performance profile is a runtime
  *distribution*: with a single solver the reference time is its own best.
  Calling it a Dolan–Moré comparison would read more into the data than it
  holds.
- **Six attempts that were built, measured and reverted.** General bound
  tightening, duplicate-column merging, a raised iteration cap, the
  bucketed search's first pass, and exact steepest-edge pricing. Each is
  recorded in `NOTICE_ALGORITHMS.md` with the measurement that killed it,
  because knowing what does not work is worth keeping.

## Verifying the claims

Nothing here asks to be taken on trust:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
ctest --test-dir build            # 10 suites
bash bench/verify_clean_room.sh   # no third-party solver in the dependency graph
python3 bench/make_report.py      # every headline number, computed not transcribed
```

`verify_clean_room.sh` is the mechanical form of the project's central
claim. It scans sources and build files, inspects what the built binaries
*actually link* with `otool`/`ldd`, and lists every CMake dependency. It
exits non-zero on any hit.
