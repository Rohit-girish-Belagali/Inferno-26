#!/usr/bin/env bash
# INFERNO — extended demo for SIH 2026 PS 26119.
#
# Designed to be screen-recorded: runs continuously for roughly 15-25
# minutes with live output throughout (every stage is unbuffered, so a
# recording shows work happening rather than a frozen terminal followed by
# a dump). Nothing here is staged or pre-computed — every number is
# produced by the engine as you watch.
set -uo pipefail
cd "$(dirname "$0")/.."
B=build

hr()   { printf '\n\033[1m════════════════════════════════════════════════════════════\033[0m\n'; }
stage(){ printf '\n\033[1;33m### %s\033[0m\n' "$1"; }
note() { printf '  \033[2m%s\033[0m\n' "$1"; }

if [ ! -x "$B/solver" ]; then
  echo "Build first: cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"
  exit 1
fi

START=$(date +%s)
elapsed(){ printf '\033[2m[t+%ds]\033[0m\n' "$(( $(date +%s) - START ))"; }

hr
printf '\033[1mINFERNO — Indigenous Optimization Solver\033[0m\n'
printf 'SIH 2026 · Problem Statement 26119 · Team INFERNO\n'
printf 'Extended demonstration run — all output live, nothing pre-computed.\n'
hr

stage "STAGE 1/9 — Correctness gate: the full test suite"
note "10 suites covering LP, QP, presolve/postsolve, linear algebra, the"
note "C ABI (compiled as C), the Python bindings and every model."
(cd $B && ctest --output-on-failure 2>&1) | tail -20
elapsed

stage "STAGE 2/9 — Sovereignty: no third-party solver, checked mechanically"
./bench/verify_clean_room.sh
elapsed

stage "STAGE 3/9 — The PS 26119 problem battery"
note "Industrial models plus degeneracy, ill-conditioning and scale stress."
$B/battery --long
elapsed

stage "STAGE 4/9 — MILP scale characterisation"
note "Ramps binary count under a hard time limit and reports what actually"
note "happened, including the sizes where optimality is NOT proved."
$B/mip_scale 30
elapsed

stage "STAGE 5/9 — Refinery crude blending, MRPL's own domain"
$B/refinery_demo
elapsed

stage "STAGE 6/9 — Dual simplex over the Netlib set"
note "A second, independent solve path. It supports only instances with a"
note "trivial dual-feasible start and reports NUMERICAL_ERROR on the rest"
note "rather than guessing — so the score is far below the primal path's"
note "by design, not by failure."
python3 -u bench/run_netlib.py --solver dual 2>&1 | tail -25
elapsed

stage "STAGE 7/9 — Presolve on/off equivalence"
note "BUILD_PLAN_V2.md makes this a standing rule: presolve must not change"
note "the answer. Same instances, presolve disabled, compared to Stage 8."
python3 -u bench/run_netlib.py --solver revised 2>&1 | tail -12
elapsed

stage "STAGE 8/9 — Full Netlib benchmark with presolve, all 93 instances, live"
note "This is the long stage. Each line is one real industrial LP, solved"
note "and then independently verified. Expect roughly 10 minutes."
python3 -u bench/run_netlib.py --solver revised --presolve
elapsed

stage "STAGE 9/9 — Benchmark report and performance profile"
python3 bench/make_report.py
elapsed

hr
printf '\033[1mRun complete.\033[0m Total elapsed: %d seconds.\n' "$(( $(date +%s) - START ))"
note "Every figure above was computed during this run. Reproduce with README.md."
hr
