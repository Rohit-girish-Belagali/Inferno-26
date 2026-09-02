# Project Status & Checkpoint Assessment

## Completed & Active Phase Summary

### Phase 1.1 — Foundation (days 1–3): ✅ Complete
- Reader, checker, throwaway dense simplex, bench harness.
- Checkpoint 1.1 met (`afiro` matches published optimum).

### Phase 1.2 — Sparse LU (days 4–10): ✅ Complete
- Markowitz LU, Gilbert-Peierls solves, scaling, refactorization policy, iterative refinement.
- One deliberate substitution: PFI eta update instead of true Forrest-Tomlin (documented, not hidden).
- Day-10 gate metrics all met (factors correctly, updates within 1e-8, ≥20x faster than refactoring).

### Phase 2.1 — Revised Simplex (days 11–25): 🟡 Partial (Short of Gate, much closer)
- Real revised simplex built — Devex pricing (with a safe deterministic tie-breaking perturbation, added this session), scaling, two-pass ratio test, and a greedy weighted-matching crash basis (also added this session).
- **80/93 Netlib instances verified with --presolve, zero checker-rejected claims** — the session's best confirmed number, up from a 77/93 starting baseline. Only 4 non-converging instances left (`pilot87`/`scsd1`/`scsd8` numerical error, `tuff` iteration limit), down from 18 at session start.
- The plan's own Checkpoint 2.1 (day 25) explicitly wants 98/98 at 1e-6 — still a real gap, meaningfully narrower than before.
- Missing: steepest-edge pricing, full Harris bound-relaxation, bound-flipping dual ratio test. Cost perturbation is done, but deliberately not the textbook form — see NOTICE_ALGORITHMS.md for why perturbing the objective itself was judged too risky to certify correct.
- Built ahead of schedule: dual simplex (partial — only instances with a trivial dual-feasible start; a general dual phase 1 isn't implemented) and a sequential solve-manager fallback.

### Phase 2.2 — Presolve (days 26–32): 🟡 Partial, wider than before
- 5 of ~8 reductions now: fixed-variable, empty-column (original two) plus redundant-row, singleton-row, and free-column-singleton removal (added this session, each with real postsolve dual/primal recovery, not just bound narrowing).
- Missing: forcing rows, dominated columns/dual fixing, duplicate row/column detection, coefficient tightening. General (non-removing) bound tightening AND duplicate-column merging were both tried and reverted — the first is sound for feasibility but breaks the checker without a harder dual-reconciliation scheme than this session took on; the second fired often on real data (1943 duplicate columns across 34 Netlib instances) and got two real bugs fixed, but a third, deeper cross-reduction ordering issue surfaced that would need unifying two separate postsolve dependency passes to fix properly — reverted rather than shipped on unverified confidence. Both documented in NOTICE_ALGORITHMS.md, not silently dropped.
- **Full Netlib (`--presolve`, revised solver, with crash basis + perturbed pricing): 80/93 verified, zero checker-rejected claims.** Fixed `greenbeb` (previously wouldn't converge in 1158s+) and `pilotnov`; net effect across all additions this session was clearly positive.
- Phase 2 gate ("30% size reduction") not formally measured yet, but the reduction set is meaningfully larger now.

### Phase 3 — GPU / PDLP (days 20–40): ⏸️ Not Started
- Nothing. Explicitly blocked (no GPU hardware on dev machine, per plan note). UI/dashboard features simulated.

### Phase 4 — MILP / QP (days 33–60): ⏸️ Not Started
- Nothing. No branch-and-bound, cuts, heuristics, or QP/ADMM.

### Phase 5 — Refinery Models, Benchmark Report, Packaging (days 55–75): 🟡 Early / Partial
- The actual refinery/crude-blending models and Dolan-Moré benchmark report aren't built.
- What is done early: the demo dashboard and bench-results report (wireframe modules 5 & 6), which overlaps with Phase 5's "demo" intent even though built out of order.

---

## Bottom Line

Solid on Phase 1 (both gates cleared honestly). Phase 2 has real, working substance but is short of its own checkpoint — the plan itself says a missed day-25 checkpoint should trigger narrowing scope to "LP + GPU path + refinery models" rather than continuing to chase MILP. This session moved Phase 2.1+2.2's combined score from 77/93 to 80/93 (zero checker-rejected claims throughout every change, verified against the full Netlib set at each step, not just spot-checked) and cut non-converging instances from 18 to 4 — real progress toward the 98/98 gate, but still short of it, and two attempted additions (general bound tightening, duplicate-column merging) were found unsound-to-ship and honestly reverted rather than pushed through. Phases 3–5 are essentially untouched. Ahead of a typical day-2 pace content-wise, but not yet past its own first real gate.
