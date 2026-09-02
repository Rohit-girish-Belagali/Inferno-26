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
| `simplex/revised_simplex.*` | Bounded-variable primal revised simplex; Devex pricing with a Bland's-rule fallback after sustained degenerate pivots; a two-pass ratio test (minimum step length, then largest-\|rate\| among ties); geometric-mean scaling (`la/scaling.*`) applied around the whole solve | Bazaraa, Jarvis & Sherali, *Linear Programming and Network Flows*, ch. 7 (revised simplex, bounded variables); Harris, "Pivot selection methods of the Devex LP code," *Math. Programming* 5, 1973 (Devex weights, and the same paper's motivation for the two-pass ratio-test structure used here — not its full bound-relaxation machinery); Bland 1977 (anti-cycling fallback) | Permanent — this is the real solve path, not throwaway. No steepest-edge, full Harris bound relaxation, dual simplex, cost perturbation, or crash basis yet — see file header. **Measured honestly across three increments** on the Netlib set: Dantzig alone 72/93; +Devex 71/93 (net neutral — fixed 4, regressed 5, big per-instance wins like d2q06c's 117786→12880 iterations, but no fix for the underlying ill-conditioning risk); +scaling (which existed since Phase 1.2 but was never actually wired into any solve path until this) 74/93 with **zero** checker-rejected "optimal" claims, versus 1 before. Scaling attacked the actual root cause (an unscaled matrix amplifying through elimination) rather than only the symptom two earlier defensive checks (`MaxBoundViolation`, before trusting an unbounded *or* optimal conclusion) were added to catch. |
| `presolve/presolve.*` | Fixed-variable substitution and empty-column removal, to a fixpoint | Standard presolve reductions — see e.g. Andersen & Andersen, "Presolving in linear programming," *Math. Programming* 71, 1995 | Partial — only the two reductions that don't remove a row (and so don't need row-dual recovery through postsolve) are implemented. Empty/singleton row removal, forcing/redundant rows, bound tightening, dominated columns, dual fixing, duplicate row/column detection, coefficient tightening are BUILD_PLAN_V2.md Phase 2.2 checklist items not yet done — see file header. **Measured on the full Netlib set** (`--presolve`, opt-in, not the CLI/bench default): 75/93 vs. 74/93 without, zero checker-rejected claims either way — net positive but not a clean win. It fixed the entire "pilot" family (pilot, pilot.ja, pilot.we, pilot4, pilotnov, maros-r7, perold — 7 instances that timed out or errored without it) but regressed degen3 and grow22 (both already numerically marginal — grow22 in particular has flipped pass/fail across several earlier changes this session too, suggesting real sensitivity in that instance rather than a presolve-specific bug). Left opt-in rather than promoted to default: the partial reduction set is explicitly incomplete, and BUILD_PLAN_V2.md's standing rule wants the full reduction set validated for on/off equivalence before it's load-bearing. |

## Third-party tools (not solver code)

`bench/tools/emps.c` is fetched verbatim from https://www.netlib.org/lp/data/emps.c
to decode the Netlib LP archive's on-disk encoding into plain MPS text. It
performs no optimization, is never linked into `inferno_core` or `solver`,
and is a data-format utility, not a solver dependency. See
`bench/tools/NOTICE.md`.

## Pending

Full Harris two-pass ratio test (bound relaxation), steepest-edge pricing,
dual simplex, the row-removing presolve reductions (empty/singleton rows,
forcing rows, dominated columns, duplicate rows, coefficient tightening),
PDLP, branch-and-bound, cuts, ADMM for QP, true Forrest-Tomlin — each gets
its citation added here when that module is written, per the core reading
list in `BUILD_PLAN_V2.md`.
