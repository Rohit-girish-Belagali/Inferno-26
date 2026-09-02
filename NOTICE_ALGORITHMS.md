# Algorithm citations

Standing rule (`BUILD_PLAN_V2.md` #3): every algorithm gets a citation here.
Implementing a published method from its description is expected; copying
code from any existing solver is not — not once, not anywhere.

| Component | Algorithm | Citation | Status |
|---|---|---|---|
| `simplex/dense_simplex.*` | Two-phase bounded-variable primal simplex, dense tableau, Bland's rule anti-cycling | Standard method — see e.g. Bazaraa, Jarvis & Sherali, *Linear Programming and Network Flows*, ch. 4 (bounded variables); Bland, "New finite pivoting rules for the simplex method," *Math. of OR* 2(2), 1977 | **Throwaway** — Phase 1.1 proof-of-life only, to be deleted once Phase 2.1's revised simplex passes Checkpoint 2.1. Not intended to survive the project. |
| `io/mps_reader.*` | MPS free-form parsing | Format only, not an algorithm — no citation needed. Netlib RANGES-section sign convention follows the description in Maros, *Computational Techniques of the Simplex Method*, App. B | Permanent |
| `checker/checker.*` | Primal/dual residual and complementarity-gap computation | Standard LP optimality (KKT) conditions — see e.g. Nocedal & Wright, *Numerical Optimization*, ch. 12 | Permanent |

## Third-party tools (not solver code)

`bench/tools/emps.c` is fetched verbatim from https://www.netlib.org/lp/data/emps.c
to decode the Netlib LP archive's on-disk encoding into plain MPS text. It
performs no optimization, is never linked into `inferno_core` or `solver`,
and is a data-format utility, not a solver dependency. See
`bench/tools/NOTICE.md`.

## Pending

Every entry above phases 2 onward (Markowitz LU, Forrest-Tomlin, Harris
ratio test, Devex/steepest-edge pricing, presolve reductions, PDLP,
branch-and-bound, cuts, ADMM for QP, ...) gets its citation added here when
that module is written, per the core reading list in `BUILD_PLAN_V2.md`.
