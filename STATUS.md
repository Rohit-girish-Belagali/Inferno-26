# Project Status & Checkpoint Assessment

## Completed & Active Phase Summary

### Phase 1.1 — Foundation (days 1–3): ✅ Complete
- Reader, checker, throwaway dense simplex, bench harness.
- Checkpoint 1.1 met (`afiro` matches published optimum).

### Phase 1.2 — Sparse LU (days 4–10): ✅ Complete
- Markowitz LU, Gilbert-Peierls solves, scaling, refactorization policy, iterative refinement.
- One deliberate substitution: PFI eta update instead of true Forrest-Tomlin (documented, not hidden).
- Day-10 gate metrics all met (factors correctly, updates within 1e-8, ≥20x faster than refactoring).

### Phase 2.1 — Revised Simplex (days 11–25): 🟡 Partial (Short of Gate)
- Real revised simplex built — Devex pricing, scaling, two-pass ratio test.
- 74–75/93 Netlib instances verified, zero false "optimal" claims.
- The plan's own Checkpoint 2.1 (day 25) explicitly wants 98/98 at 1e-6 — this is a real gap, not close.
- Missing: steepest-edge pricing, full Harris bound-relaxation, cost perturbation, crash basis — all named as the reason the remaining ~18 hard instances (`pilot87`, `dfl001`, `degen3`-class problems) don't solve.
- Built ahead of schedule: dual simplex (partial — only instances with a trivial dual-feasible start; a general dual phase 1 isn't implemented) and a sequential solve-manager fallback.

### Phase 2.2 — Presolve (days 26–32): 🟡 Early / Partial
- Only 2 of ~8 reductions (fixed-variable, empty-column) — the ones that don't need row-dual recovery.
- Missing: empty/singleton row removal, forcing rows, bound tightening, dominated columns, duplicate detection, coefficient tightening.
- Phase 2 gate ("30% size reduction") isn't met with just these two.

### Phase 3 — GPU / PDLP (days 20–40): ⏸️ Not Started
- Nothing. Explicitly blocked (no GPU hardware on dev machine, per plan note). UI/dashboard features simulated.

### Phase 4 — MILP / QP (days 33–60): ⏸️ Not Started
- Nothing. No branch-and-bound, cuts, heuristics, or QP/ADMM.

### Phase 5 — Refinery Models, Benchmark Report, Packaging (days 55–75): 🟡 Early / Partial
- The actual refinery/crude-blending models and Dolan-Moré benchmark report aren't built.
- What is done early: the demo dashboard and bench-results report (wireframe modules 5 & 6), which overlaps with Phase 5's "demo" intent even though built out of order.

---

## Bottom Line

Solid on Phase 1 (both gates cleared honestly). Phase 2 has real, working substance but is short of its own checkpoint — the plan itself says a missed day-25 checkpoint should trigger narrowing scope to "LP + GPU path + refinery models" rather than continuing to chase MILP. Phases 3–5 are essentially untouched. Ahead of a typical day-2 pace content-wise, but not yet past its own first real gate.
