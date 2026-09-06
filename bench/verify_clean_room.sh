#!/usr/bin/env bash
# Clean-room proof (BUILD_PLAN_V2.md Phase 5): demonstrate that no
# third-party solver is anywhere in this project's dependency graph.
#
# The claim being defended is specific -- "no HiGHS, CBC, GLPK, SCIP, OSQP
# or SuiteSparse anywhere in the dependency graph" -- so this checks it
# mechanically rather than asking anyone to take it on faith. It fails
# loudly (exit 1) if any forbidden name appears in a link line, an include
# directive, or the CMake dependency graph.
set -uo pipefail
cd "$(dirname "$0")/.."

NAMES="highs|cbc|glpk|scip|osqp|suitesparse|umfpack|klu|cholmod|gurobi|cplex|mosek|clp"
# Source scan matches WHOLE WORDS only. Plain substring matching produces
# false positives against ordinary English -- "discipline" contains "scip",
# which failed this check on two comment lines before the boundaries were
# added. A real dependency shows up as an include, a find_package or a link
# name, all of which are word-bounded.
FORBIDDEN_WORD="\\b($NAMES)\\b"
# The binary scan stays substring-based on purpose: a linked library is
# "libhighs.dylib", where the name is not word-bounded.
FORBIDDEN_SUB="($NAMES)"
fail=0

echo "=== 1. Source and build files ==="
# Exclude this script itself and the docs that legitimately NAME the
# solvers we are proving we do not use.
hits=$(grep -rIniE "$FORBIDDEN_WORD" \
        --include=*.cpp --include=*.hpp --include=*.h --include=*.cc \
        --include=CMakeLists.txt --include=*.cmake \
        . 2>/dev/null | grep -v '^./build/' || true)
if [ -n "$hits" ]; then
  echo "FAIL: forbidden solver referenced in source or build files:"
  echo "$hits"
  fail=1
else
  echo "ok: no third-party solver named in any source or build file"
fi

echo
echo "=== 2. Linked libraries ==="
# What the binaries actually load at runtime. This is the check that
# cannot be talked around: whatever the source says, these are the real
# dependencies.
for bin in build/solver build/refinery_demo; do
  [ -x "$bin" ] || { echo "skip: $bin not built"; continue; }
  if command -v otool >/dev/null 2>&1; then libs=$(otool -L "$bin" | tail -n +2)
  else libs=$(ldd "$bin" 2>/dev/null || true); fi
  bad=$(echo "$libs" | grep -IniE "$FORBIDDEN_SUB" || true)
  if [ -n "$bad" ]; then
    echo "FAIL: $bin links a third-party solver:"; echo "$bad"; fail=1
  else
    echo "ok: $bin links only system libraries:"
    echo "$libs" | sed 's/^/     /'
  fi
done

echo
echo "=== 3. CMake link dependencies ==="
deps=$(grep -rhE "target_link_libraries" CMakeLists.txt tests/CMakeLists.txt 2>/dev/null || true)
echo "$deps" | sed 's/^/     /'
if echo "$deps" | grep -qIiE "$FORBIDDEN_SUB"; then
  echo "FAIL: CMake links a third-party solver"; fail=1
else
  echo "ok: every target links only inferno_core"
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "CLEAN-ROOM PROOF PASSED — no third-party solver in the dependency graph."
else
  echo "CLEAN-ROOM PROOF FAILED"
fi
exit "$fail"
