#!/usr/bin/env bash
# INFERNO — end-to-end demo for SIH 2026 Problem Statement 26119
# ("Indigenous GPU-Accelerated Optimization Solver").
#
# Walks the analyst user flow exactly as the proposal deck describes it:
#
#   Model in (MPS) -> Parse & Validate -> Presolve -> Scale
#     -> Concurrent Solve Manager -> Exact Postsolve
#     -> Independent Solution Checker -> Report
#
# and then shows the same engine driven from each of the three modeling
# APIs the PS names (MPS, C++/C ABI, Python).
#
# It finishes with a capability scorecard against the PS's six modules,
# including the ones that are NOT built. A demo that hides its gaps is
# worth less than one that states them, because the first question from a
# panel is always the gap.
set -uo pipefail
cd "$(dirname "$0")/.."
B=build

hr()   { printf '\n\033[1m%s\033[0m\n' "────────────────────────────────────────────────────────────"; }
step() { printf '\n\033[1;36m▶ %s\033[0m\n' "$1"; }
note() { printf '  \033[2m%s\033[0m\n' "$1"; }

if [ ! -x "$B/solver" ]; then
  echo "Build first:  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"
  exit 1
fi

hr
printf '\033[1mINFERNO — Indigenous Optimization Solver\033[0m\n'
printf 'SIH 2026 · Problem Statement 26119 · Team INFERNO\n'
note "Engine version $($B/solver --version 2>/dev/null || echo 0.5.0) · no third-party solver linked"

# ---------------------------------------------------------------- 1
step "1. MODEL IN — parse and validate a real industrial LP (MPS)"
MPS=bench/netlib/mps/afiro.mps
if [ ! -f "$MPS" ]; then
  note "Netlib set not present; fetch it with bench/download_netlib.sh"
  note "Skipping the MPS leg of the demo."
else
  note "Reading $MPS"
  $B/solver "$MPS" --solver revised | sed 's/^/  /'
  note "The 'checker:' line is an INDEPENDENT verification — it recomputes"
  note "the residuals from the raw problem and shares no code with the solver."
fi

# ---------------------------------------------------------------- 2
step "2. PRESOLVE + POSTSOLVE — shrink, solve, then reverse every reduction"
if [ -f "$MPS" ]; then
  $B/solver "$MPS" --solver revised --presolve | sed 's/^/  /'
  note "Six reduction kinds run to a fixpoint; postsolve reverses each one"
  note "exactly, recovering duals for rows presolve removed."
fi

# ---------------------------------------------------------------- 3
step "3. CONCURRENT SOLVE MANAGER — race solve paths, report the winner"
if [ -f "$MPS" ]; then
  $B/solver "$MPS" --solver managed | sed 's/^/  /'
  note "Reports which path produced the answer. Honest scope: this is a"
  note "sequential fallback chain, not true concurrency — see solve_manager.hpp."
fi

# ---------------------------------------------------------------- 4
step "4. THE INDUSTRY MODEL — refinery crude blending (MRPL's own domain)"
$B/refinery_demo | sed 's/^/  /'

# ---------------------------------------------------------------- 5
step "5. BREADTH — the same engine on other sectors"
note "Transportation (supply chain) and power-system economic dispatch"
$B/tests/breadth_test 2>/dev/null | tail -3 | sed 's/^/  /'

# ---------------------------------------------------------------- 6
step "6. QUADRATIC PROGRAMMING — convex QP via ADMM"
$B/tests/qp_test 2>/dev/null | tail -4 | sed 's/^/  /'

# ---------------------------------------------------------------- 7
step "7. MODELING APIs — the PS names MPS, C++ and Python"
note "C ABI (a C consumer, compiled as C):"
$B/tests/api_test 2>/dev/null | tail -3 | sed 's/^/    /'
note "Python:"
python3 - <<'PY' 2>&1 | sed 's/^/    /'
import sys, os
sys.path.insert(0, "python")
import inferno
p = inferno.Problem.from_mps("bench/netlib/mps/afiro.mps") if os.path.exists("bench/netlib/mps/afiro.mps") else None
if p is None:
    print("Netlib set not present; skipping.")
else:
    r = p.solve(presolve=True)
    print(f"import inferno -> {r.status_name}  objective={r.objective:.6f}  checker_passed={r.checker_passed}")
    p.close()
PY

# ---------------------------------------------------------------- 8
step "8. BENCHMARKING & ANALYTICS — the whole Netlib set"
python3 bench/make_report.py 2>/dev/null | sed -n '1,14p' | sed 's/^/  /'

# ---------------------------------------------------------------- 9
step "9. SOVEREIGNTY — the central claim, checked mechanically"
./bench/verify_clean_room.sh 2>/dev/null | grep -E "^(ok|FAIL|CLEAN-ROOM)" | sed 's/^/  /'

# ---------------------------------------------------------------- 10
hr
printf '\033[1mCAPABILITY SCORECARD vs the six PS 26119 modules\033[0m\n\n'
cat <<'TABLE'
  1. Solver Core Engine (LP/MILP/QP)   PARTIAL  LP 91/93 Netlib, 0 wrong answers.
                                                QP via ADMM. MILP NOT built —
                                                dropped by the plan's own day-25
                                                rule, which says a correct LP
                                                solver beats a broken MILP one.

  2. GPU Acceleration Layer (CUDA)     NOT BUILT  No NVIDIA GPU on this machine,
                                                so no CUDA port and NO GPU
                                                SPEEDUP IS CLAIMED. The CPU PDLP
                                                reference the port would mirror
                                                IS built and solves 4/15 small
                                                instances to full accuracy.

  3. AI-Assisted Presolve & Heuristics PARTIAL  Presolve is real: 6 reductions
                                                with exact postsolve. The
                                                "AI-assisted" part is NOT built,
                                                and is not simulated.

  4. Modeling APIs (Python/C++/MPS)    DONE     All three, all exercised above
                                                and all covered by tests.

  5. Benchmarking & Analytics          DONE     Full Netlib harness scored by an
                                                independent checker, plus a
                                                report and performance profile.

  6. Dashboard & Cloud Deploy          PARTIAL  Dashboard pages exist on real
                                                measured data. Cloud deploy NOT
                                                built.
TABLE
printf '\n'
note "Everything above is reproducible: see README.md."
hr
