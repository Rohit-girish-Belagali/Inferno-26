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
- **91/93 Netlib instances verified with --presolve, zero checker-rejected claims, zero numerical errors**, all 91 solving in 61.3s total — up from a 77/93 starting baseline. Only `dfl001` and `greenbea` remain, both `TIMEOUT`.
- **The `NUMERICAL_ERROR` category is empty and no instance reports a wrong answer.** The two remaining failures are pure wall-clock timeouts on the largest instances (`dfl001`, `greenbea`), not correctness or robustness failures.
- The plan's own Checkpoint 2.1 (day 25) explicitly wants 98/98 at 1e-6 — still a real gap, but the nature of it has changed: what remains is a speed and degeneracy problem, not a robustness one.
- **Refactorization interval now scales with basis size** (was a flat 100 eta updates): component timing showed refactorization was 77% of pilot87's solve. Scaling it newly solved `degen3`, `fit2p` and `pilot87` — three of the hard instances, two of them named in the plan itself — with zero regressions, taking the score 81/93 → 84/93. Later retuned to m/20 alongside the bucketed LU pivot search, which only work as a pair (see NOTICE_ALGORITHMS.md): together they cut the whole verified set from ~93s to 63.5s, with `pilot87` 56.9s → 14.9s, `degen3` 45.6s → 12.0s, `fit2p` 46.3s → 7.6s.
- **Performance: 1.57x faster** across the 81 verified instances (163.6s → 104.2s total), from sparse-scattering the Devex weight update — profiling put ~54% of a hard solve in that one loop, which was heap-allocating a vector per nonbasic column per pivot. Biggest wins `fit2d` 23.3s → 3.8s, `scsd8` 2.8s → 0.8s. Same objectives and same iteration counts, so this was purely per-iteration cost.
- **Where the remaining gap actually is:** the degeneracy problem is solved — the Harris ratio test with adaptive tolerance expansion took the score 84/93 → 91/93, clearing every stalling instance (`tuff`, `wood1p`, `modszk1`, `cycle`, `d6cube`). What is left is raw throughput on the two largest models. Exact steepest-edge pricing was tried for those and reverted — it regressed the hard instances badly because the exact recurrence needs exact initial weights, which cost one FTRAN per column (see NOTICE_ALGORITHMS.md). Making it work needs a reference-framework reset scheme, not just the recurrence.
- Missing: steepest-edge pricing (attempted, reverted — see NOTICE_ALGORITHMS.md), bound-flipping dual ratio test. The full Harris ratio test with tolerance expansion is now done. Cost perturbation is done, but deliberately not the textbook form — see NOTICE_ALGORITHMS.md for why perturbing the objective itself was judged too risky to certify correct.
- Built ahead of schedule: dual simplex (partial — only instances with a trivial dual-feasible start; a general dual phase 1 isn't implemented) and a sequential solve-manager fallback.

### Phase 2.2 — Presolve (days 26–32): 🟡 Partial, wider than before
- 5 of ~8 reductions now: fixed-variable, empty-column (original two) plus redundant-row, singleton-row, and free-column-singleton removal (added this session, each with real postsolve dual/primal recovery, not just bound narrowing).
- Missing: forcing rows, dominated columns/dual fixing, duplicate row/column detection, coefficient tightening. General (non-removing) bound tightening AND duplicate-column merging were both tried and reverted — the first is sound for feasibility but breaks the checker without a harder dual-reconciliation scheme than this session took on; the second fired often on real data (1943 duplicate columns across 34 Netlib instances) and got two real bugs fixed, but a third, deeper cross-reduction ordering issue surfaced that would need unifying two separate postsolve dependency passes to fix properly — reverted rather than shipped on unverified confidence. Both documented in NOTICE_ALGORITHMS.md, not silently dropped.
- **Full Netlib (`--presolve`, revised solver): 91/93 verified, zero checker-rejected claims, zero numerical errors.** Fixed `greenbeb` (previously wouldn't converge in 1158s+) and `pilotnov`; net effect across all additions this session was clearly positive.
- Phase 2 gate ("30% size reduction") not formally measured yet, but the reduction set is meaningfully larger now.

### Phase 3 — GPU / PDLP (days 20–40): ⏸️ Not Started
- Nothing. Explicitly blocked (no GPU hardware on dev machine, per plan note). UI/dashboard features simulated.

### Phase 4 — MILP / QP (days 33–60): ⏸️ Not Started
- Nothing. No branch-and-bound, cuts, heuristics, or QP/ADMM.

### Phase 5 — Refinery Models, Benchmark Report, Packaging (days 55–75): 🟡 Early / Partial
- The actual refinery/crude-blending models and Dolan-Moré benchmark report aren't built.
- What is done early: the demo dashboard and bench-results report (wireframe modules 5 & 6), which overlaps with Phase 5's "demo" intent even though built out of order.

---

## Known issue for whoever picks this up

`la_test` is intermittently red, and the cause is understood: it asserts a
**wall-clock ratio** (basis update at least 20x faster than a fresh
refactorization) inside a unit test. Measured on an idle machine it lands at
26–32x, so the margin is real — but running it while a full Netlib bench is
going pushes it under the 20x floor and the suite fails. Every observed
failure was that test, under exactly those conditions, and it passes on
re-run. It is a test-design problem (a performance assertion in a
correctness suite), not a solver regression: re-run `ctest` on a quiet
machine before believing a red result, and consider making that assertion
either load-tolerant or a separate benchmark target.

---

## Bottom Line

Solid on Phase 1 (both gates cleared honestly). Phase 2 has real, working substance but is short of its own checkpoint — the plan itself says a missed day-25 checkpoint should trigger narrowing scope to "LP + GPU path + refinery models" rather than continuing to chase MILP. This session moved Phase 2.1+2.2's combined score from 77/93 to 91/93 (zero checker-rejected claims throughout every change, verified against the full Netlib set at each step, not just spot-checked) and — more significant than the count — eliminated the `NUMERICAL_ERROR` category entirely, so the solver no longer abandons any instance with an error. Everything still unsolved is a speed/iteration-budget problem, which is a far better position to attack the 98/98 gate from than a robustness one. Two attempted additions (general bound tightening, duplicate-column merging) were found unsound-to-ship and honestly reverted rather than pushed through. Phases 3–5 are essentially untouched. Ahead of a typical day-2 pace content-wise, but not yet past its own first real gate.
