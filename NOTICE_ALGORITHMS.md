# Algorithm citations

Standing rule (`BUILD_PLAN_V2.md` #3): every algorithm gets a citation here.
Implementing a published method from its description is expected; copying
code from any existing solver is not — not once, not anywhere.

| Component | Algorithm | Citation | Status |
|---|---|---|---|
| `simplex/dense_simplex.*` | Two-phase bounded-variable primal simplex, dense tableau, Bland's rule anti-cycling | Standard method — see e.g. Bazaraa, Jarvis & Sherali, *Linear Programming and Network Flows*, ch. 4 (bounded variables); Bland, "New finite pivoting rules for the simplex method," *Math. of OR* 2(2), 1977 | **Throwaway** — Phase 1.1 proof-of-life only, to be deleted once Phase 2.1's revised simplex passes Checkpoint 2.1. Not intended to survive the project. |
| `io/mps_reader.*` | MPS free-form parsing, with a strict fixed-column fallback | Format only, not an algorithm — no citation needed. Netlib RANGES-section sign convention and the fixed-column field layout (columns 2-3/5-12/15-22/25-36/40-47/50-61) follow the description in Maros, *Computational Techniques of the Simplex Method*, App. B | Permanent |
| `checker/checker.*` | Primal/dual residual and complementarity-gap computation | Standard LP optimality (KKT) conditions — see e.g. Nocedal & Wright, *Numerical Optimization*, ch. 12 | Permanent |
| `la/scaling.*` | Geometric-mean scaling + infinity-norm equilibration | Suhl & Suhl 1990, "A fast LU update for linear programming" | Permanent |
| `la/markowitz_lu.*` | Markowitz-count pivot selection with threshold partial pivoting | Markowitz, "The elimination form of the inverse and its application to linear programming," *Management Science* 3(3), 1957; threshold pivoting and elimination bookkeeping per Suhl & Suhl 1990 | Permanent, MVP pivot search (see file header) |
| `la/lu_solve.*` | Gilbert-Peierls sparse triangular solve (DFS reachability) | Gilbert & Peierls, "Sparse partial pivoting in time proportional to arithmetic operations," *SIAM J. Sci. Stat. Comput.* 9(5), 1988 | Permanent |
| `la/basis_factorization.*` | Product-form-of-the-inverse (PFI) eta update, **not** Forrest-Tomlin | Classical technique — see e.g. Maros, *Computational Techniques of the Simplex Method*, ch. 3. **This is a deliberate, documented substitution**: BUILD_PLAN_V2.md's Phase 1.2 checklist names Forrest & Tomlin 1972's LU-form bump/Hessenberg update specifically. That algorithm's exact permutation-and-elimination procedure was judged too easy to get subtly wrong from memory without a reference to check against, so the simpler, fully-rederived PFI eta update was implemented instead, paired with a refactorization policy to bound the eta-chain growth PFI alone doesn't guarantee against the way true Forrest-Tomlin does. See the file header comment. | Interim — upgrading to true Forrest-Tomlin is tracked, not done |

## Third-party tools (not solver code)

`bench/tools/emps.c` is fetched verbatim from https://www.netlib.org/lp/data/emps.c
to decode the Netlib LP archive's on-disk encoding into plain MPS text. It
performs no optimization, is never linked into `inferno_core` or `solver`,
and is a data-format utility, not a solver dependency. See
`bench/tools/NOTICE.md`.

## Pending

Every entry above Phase 2 onward (revised simplex, Harris ratio test,
Devex/steepest-edge pricing, presolve reductions, PDLP, branch-and-bound,
cuts, ADMM for QP, true Forrest-Tomlin, ...) gets its citation added here
when that module is written, per the core reading list in
`BUILD_PLAN_V2.md`.
