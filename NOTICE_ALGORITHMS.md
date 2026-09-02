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
| `simplex/revised_simplex.*` | Bounded-variable primal revised simplex; Devex pricing with a Bland's-rule fallback after sustained degenerate pivots; a two-pass ratio test (minimum step length, then largest-\|rate\| among ties) | Bazaraa, Jarvis & Sherali, *Linear Programming and Network Flows*, ch. 7 (revised simplex, bounded variables); Harris, "Pivot selection methods of the Devex LP code," *Math. Programming* 5, 1973 (Devex weights, and the same paper's motivation for the two-pass ratio-test structure used here — not its full bound-relaxation machinery); Bland 1977 (anti-cycling fallback) | Permanent — this is the real solve path, not throwaway. No steepest-edge, full Harris bound relaxation, dual simplex, cost perturbation, or crash basis yet — see file header. **Devex vs. Dantzig, measured honestly**: replacing Dantzig with Devex was net roughly neutral on the Netlib set (72/93 → 71/93 verified) — it fixed 4 instances (including a real wrong-answer bug on scsd6) and regressed 5 others (bore3d, fit1p, grow15, maros-r7, pilot4: previously-fine instances now hitting NUMERICAL_ERROR or timing out). The wins are real (d2q06c: 117786 iterations → 12880, a 9x reduction) but Devex alone doesn't fix the underlying ill-conditioning risk from having no Harris bound-relaxation or anti-degeneracy handling yet — that's still the more direct fix for the regressions, tracked as the next increment. |

## Third-party tools (not solver code)

`bench/tools/emps.c` is fetched verbatim from https://www.netlib.org/lp/data/emps.c
to decode the Netlib LP archive's on-disk encoding into plain MPS text. It
performs no optimization, is never linked into `inferno_core` or `solver`,
and is a data-format utility, not a solver dependency. See
`bench/tools/NOTICE.md`.

## Pending

Full Harris two-pass ratio test (bound relaxation), Devex/steepest-edge
pricing, dual simplex, presolve reductions, PDLP, branch-and-bound, cuts,
ADMM for QP, true Forrest-Tomlin — each gets its citation added here when
that module is written, per the core reading list in `BUILD_PLAN_V2.md`.
