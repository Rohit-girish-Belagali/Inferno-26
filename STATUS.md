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
- Real revised simplex built — Devex pricing (with a safe deterministic tie-breaking perturbation), scaling, two-pass ratio test, a greedy weighted-matching crash basis, and singular-basis repair with phase restart (all added this session).
- **81/93 Netlib instances verified with --presolve, zero checker-rejected claims** — up from a 77/93 starting baseline.
- **The `NUMERICAL_ERROR` category is now empty.** All 12 remaining non-verified instances are time/iteration budget issues (11 `TIMEOUT`, 1 `ITERATION_LIMIT`), not correctness failures — the solver no longer gives up on any instance with an error. `pilot87` in particular now solves correctly to its published optimum (301.71, checker PASS) but needs ~110s against the bench's 60s per-instance cap, so it still scores as `TIMEOUT`.
- The plan's own Checkpoint 2.1 (day 25) explicitly wants 98/98 at 1e-6 — still a real gap, but the nature of it has changed: what remains is a SPEED problem (and `tuff`'s iteration limit), not a robustness one.
- **Performance: 1.57x faster** across the 81 verified instances (163.6s → 104.2s total), from sparse-scattering the Devex weight update — profiling put ~54% of a hard solve in that one loop, which was heap-allocating a vector per nonbasic column per pivot. Biggest wins `fit2d` 23.3s → 3.8s, `scsd8` 2.8s → 0.8s. Same objectives and same iteration counts, so this was purely per-iteration cost.
- **Where the remaining gap actually is, measured not guessed:** after that fix, profiling shows the cost is now the pricing loop, and the binding constraint is ITERATION COUNT, not per-iteration cost. `pilot87` solves correctly (301.71, checker PASS) in 97s but needs 22323 iterations. The same sparse-scatter trick was considered for pricing and rejected on measurement: `rho` (a BTRAN of a unit vector) is genuinely sparse, but `y` (a BTRAN of the cost vector) measured 55–69% dense and rising, so it would have bought little while perturbing pivot paths. **Steepest-edge pricing is therefore the clear next lever** — it attacks iteration count directly, which is what these instances actually need.
- Missing: steepest-edge pricing, full Harris bound-relaxation, bound-flipping dual ratio test. Cost perturbation is done, but deliberately not the textbook form — see NOTICE_ALGORITHMS.md for why perturbing the objective itself was judged too risky to certify correct.
- Built ahead of schedule: dual simplex (partial — only instances with a trivial dual-feasible start; a general dual phase 1 isn't implemented) and a sequential solve-manager fallback.

### Phase 2.2 — Presolve (days 26–32): 🟡 Partial, wider than before
- 5 of ~8 reductions now: fixed-variable, empty-column (original two) plus redundant-row, singleton-row, and free-column-singleton removal (added this session, each with real postsolve dual/primal recovery, not just bound narrowing).
- Missing: forcing rows, dominated columns/dual fixing, duplicate row/column detection, coefficient tightening. General (non-removing) bound tightening AND duplicate-column merging were both tried and reverted — the first is sound for feasibility but breaks the checker without a harder dual-reconciliation scheme than this session took on; the second fired often on real data (1943 duplicate columns across 34 Netlib instances) and got two real bugs fixed, but a third, deeper cross-reduction ordering issue surfaced that would need unifying two separate postsolve dependency passes to fix properly — reverted rather than shipped on unverified confidence. Both documented in NOTICE_ALGORITHMS.md, not silently dropped.
- **Full Netlib (`--presolve`, revised solver, with crash basis + perturbed pricing + singular-basis repair): 81/93 verified, zero checker-rejected claims, zero numerical errors.** Fixed `greenbeb` (previously wouldn't converge in 1158s+) and `pilotnov`; net effect across all additions this session was clearly positive.
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

Solid on Phase 1 (both gates cleared honestly). Phase 2 has real, working substance but is short of its own checkpoint — the plan itself says a missed day-25 checkpoint should trigger narrowing scope to "LP + GPU path + refinery models" rather than continuing to chase MILP. This session moved Phase 2.1+2.2's combined score from 77/93 to 81/93 (zero checker-rejected claims throughout every change, verified against the full Netlib set at each step, not just spot-checked) and — more significant than the count — eliminated the `NUMERICAL_ERROR` category entirely, so the solver no longer abandons any instance with an error. Everything still unsolved is a speed/iteration-budget problem, which is a far better position to attack the 98/98 gate from than a robustness one. Two attempted additions (general bound tightening, duplicate-column merging) were found unsound-to-ship and honestly reverted rather than pushed through. Phases 3–5 are essentially untouched. Ahead of a typical day-2 pace content-wise, but not yet past its own first real gate.
