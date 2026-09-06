#!/usr/bin/env python3
"""Tests for the Python bindings. Run: python3 python/test_bindings.py

Registered with ctest, so a break in the bindings fails the same suite as
a break in the solver -- bindings that are not run are bindings that rot.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import inferno  # noqa: E402

failures = 0


def expect(cond, msg):
    global failures
    if cond:
        print(f"ok: {msg}")
    else:
        print(f"FAIL: {msg}", file=sys.stderr)
        failures += 1


def near(a, b, tol=1e-6):
    return abs(a - b) <= tol


expect(isinstance(inferno.version(), str), "version() returns a string")

# max 3x + 2y  ->  min -3x - 2y, s.t. x + y <= 4, x + 3y <= 6, 0 <= x,y <= 3.
# Optimum x=3, y=1, objective -11 -- the same model the C test uses, so a
# discrepancy between the two points at the bindings rather than the solver.
with inferno.Problem.from_arrays(
    num_rows=2, num_cols=2,
    col_ptr=[0, 2, 4], row_idx=[0, 1, 0, 1], values=[1, 1, 1, 3],
    obj=[-3.0, -2.0], col_lo=[0, 0], col_hi=[3, 3],
    row_lo=[-inferno.INFINITY, -inferno.INFINITY], row_hi=[4, 6],
) as p:
    expect(p.num_rows == 2 and p.num_cols == 2, "dimensions round-trip through from_arrays")
    r = p.solve()
    expect(r.status == inferno.SolveStatus.OPTIMAL, "in-memory problem solves to optimal")
    expect(r.status_name == "OPTIMAL", "status_name is readable")
    expect(near(r.objective, -11.0), f"objective is -11 (got {r.objective})")
    expect(near(r.x[0], 3.0) and near(r.x[1], 1.0), f"solution is x=3,y=1 (got {r.x})")
    expect(r.checker_passed, "independent checker accepted it")
    expect(len(r.y) == 2, "duals come back with one entry per row")

# The MPS path, against a known published optimum.
mps = "bench/netlib/mps/afiro.mps"
if os.path.exists(mps):
    with inferno.Problem.from_mps(mps) as p:
        r = p.solve(presolve=True)
        expect(r.status == inferno.SolveStatus.OPTIMAL, "afiro solves from MPS")
        expect(near(r.objective, -464.7531428571, 1e-6),
               f"afiro matches its published optimum (got {r.objective})")
        expect(r.checker_passed, "afiro passes the independent checker through Python")
else:
    print("skip: no local Netlib set")

# Errors must surface as exceptions, not silent wrong answers.
try:
    inferno.Problem.from_mps("/definitely/not/a/file.mps")
    expect(False, "a missing file raises")
except inferno.InfernoError:
    expect(True, "a missing file raises InfernoError")

# close() must be idempotent -- __del__ runs after an explicit close.
p = inferno.Problem()
p.close()
p.close()
expect(True, "close() is idempotent")

if failures:
    print(f"{failures} test(s) failed", file=sys.stderr)
    sys.exit(1)
print("all Python binding tests passed")
