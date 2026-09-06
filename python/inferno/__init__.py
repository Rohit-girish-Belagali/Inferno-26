"""Python bindings for the Inferno LP solver.

Built on the C ABI via ctypes rather than a compiled extension module,
which is a deliberate trade: ctypes needs no build step, no compiler on
the user's machine, and no matching of Python's ABI, so `import inferno`
works against a prebuilt shared library on any Python 3. A compiled
extension would be faster per call, but every call here already wraps a
solve that costs milliseconds at minimum, so the crossing overhead is
irrelevant.

    import inferno
    p = inferno.Problem.from_mps("bench/netlib/mps/afiro.mps")
    r = p.solve(presolve=True)
    print(r.status, r.objective, r.checker_passed)
"""

import ctypes
import os
from ctypes import POINTER, c_char_p, c_double, c_int, c_void_p

__all__ = ["Problem", "Result", "SolveStatus", "InfernoError", "load_library", "version"]

INFINITY = 1.0e308


class SolveStatus:
    OPTIMAL = 0
    INFEASIBLE = 1
    UNBOUNDED = 2
    ITERATION_LIMIT = 3
    NUMERICAL_ERROR = 4

    _NAMES = {
        0: "OPTIMAL",
        1: "INFEASIBLE",
        2: "UNBOUNDED",
        3: "ITERATION_LIMIT",
        4: "NUMERICAL_ERROR",
    }

    @classmethod
    def name(cls, value):
        return cls._NAMES.get(value, f"UNKNOWN({value})")


class InfernoError(RuntimeError):
    """Raised when the C ABI returns a non-OK status."""


_lib = None


def load_library(path=None):
    """Loads libinferno. Searches the usual build locations if no path is
    given, so the common case needs no configuration."""
    global _lib
    if _lib is not None and path is None:
        return _lib

    candidates = []
    if path:
        candidates.append(path)
    else:
        here = os.path.dirname(os.path.abspath(__file__))
        root = os.path.dirname(os.path.dirname(here))
        for name in ("libinferno.dylib", "libinferno.so", "inferno.dll"):
            candidates.append(os.path.join(root, "build", name))
            candidates.append(os.path.join(root, name))
        env = os.environ.get("INFERNO_LIBRARY")
        if env:
            candidates.insert(0, env)

    for cand in candidates:
        if os.path.exists(cand):
            _lib = ctypes.CDLL(cand)
            _bind(_lib)
            return _lib
    raise InfernoError(
        "could not find libinferno. Build it with "
        "`cmake --build build` or set INFERNO_LIBRARY to its path. Tried: "
        + ", ".join(candidates)
    )


def _bind(lib):
    """Declares every signature. Without this ctypes assumes int returns
    and int arguments, which silently truncates pointers on 64-bit — the
    classic way ctypes bindings crash or corrupt memory."""
    lib.inferno_problem_create.restype = c_void_p
    lib.inferno_problem_create.argtypes = []
    lib.inferno_problem_destroy.restype = None
    lib.inferno_problem_destroy.argtypes = [c_void_p]
    lib.inferno_problem_load_mps.restype = c_int
    lib.inferno_problem_load_mps.argtypes = [c_void_p, c_char_p]
    lib.inferno_problem_set.restype = c_int
    lib.inferno_problem_set.argtypes = [
        c_void_p, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_double),
        POINTER(c_double), POINTER(c_double), POINTER(c_double),
        POINTER(c_double), POINTER(c_double),
    ]
    lib.inferno_solve.restype = c_int
    lib.inferno_solve.argtypes = [c_void_p, c_int]
    for fn in ("inferno_get_solve_status", "inferno_get_iterations",
               "inferno_get_checker_passed", "inferno_get_num_rows", "inferno_get_num_cols"):
        getattr(lib, fn).restype = c_int
        getattr(lib, fn).argtypes = [c_void_p, POINTER(c_int)]
    lib.inferno_get_objective.restype = c_int
    lib.inferno_get_objective.argtypes = [c_void_p, POINTER(c_double)]
    for fn in ("inferno_get_solution", "inferno_get_duals"):
        getattr(lib, fn).restype = c_int
        getattr(lib, fn).argtypes = [c_void_p, POINTER(c_double), c_int]
    lib.inferno_status_string.restype = c_char_p
    lib.inferno_status_string.argtypes = [c_int]
    lib.inferno_version.restype = c_char_p
    lib.inferno_version.argtypes = []


def _check(lib, status):
    if status != 0:
        raise InfernoError(lib.inferno_status_string(status).decode())


def version():
    lib = load_library()
    return lib.inferno_version().decode()


class Result:
    """What a solve produced. `checker_passed` is the one to trust: it is
    the verdict of the independent checker, which recomputes the residuals
    from the raw problem rather than believing the solver."""

    def __init__(self, status, objective, iterations, checker_passed, x, y):
        self.status = status
        self.objective = objective
        self.iterations = iterations
        self.checker_passed = checker_passed
        self.x = x
        self.y = y

    @property
    def status_name(self):
        return SolveStatus.name(self.status)

    def __repr__(self):
        return (f"Result(status={self.status_name}, objective={self.objective!r}, "
                f"iterations={self.iterations}, checker_passed={self.checker_passed})")


class Problem:
    """An LP. Owns a handle to the C side and frees it deterministically —
    use it as a context manager, or let __del__ handle it."""

    def __init__(self):
        self._lib = load_library()
        self._handle = self._lib.inferno_problem_create()
        if not self._handle:
            raise InfernoError("could not allocate problem")

    @classmethod
    def from_mps(cls, path):
        p = cls()
        _check(p._lib, p._lib.inferno_problem_load_mps(p._handle, str(path).encode()))
        return p

    @classmethod
    def from_arrays(cls, num_rows, num_cols, col_ptr, row_idx, values, obj,
                    col_lo, col_hi, row_lo, row_hi):
        """Builds a problem from CSC arrays. Any sequence of numbers works;
        they are converted to C arrays here, so callers need not depend on
        numpy."""
        p = cls()

        def ints(seq):
            return (c_int * len(seq))(*[int(v) for v in seq])

        def dbls(seq):
            return (c_double * len(seq))(*[float(v) for v in seq])

        _check(p._lib, p._lib.inferno_problem_set(
            p._handle, int(num_rows), int(num_cols),
            ints(col_ptr), ints(row_idx) if len(row_idx) else ints([0]),
            dbls(values) if len(values) else dbls([0.0]),
            dbls(obj), dbls(col_lo), dbls(col_hi), dbls(row_lo), dbls(row_hi)))
        return p

    @property
    def num_rows(self):
        v = c_int()
        _check(self._lib, self._lib.inferno_get_num_rows(self._handle, ctypes.byref(v)))
        return v.value

    @property
    def num_cols(self):
        v = c_int()
        _check(self._lib, self._lib.inferno_get_num_cols(self._handle, ctypes.byref(v)))
        return v.value

    def solve(self, presolve=False):
        _check(self._lib, self._lib.inferno_solve(self._handle, 1 if presolve else 0))

        status, iters, passed = c_int(), c_int(), c_int()
        objective = c_double()
        _check(self._lib, self._lib.inferno_get_solve_status(self._handle, ctypes.byref(status)))
        _check(self._lib, self._lib.inferno_get_objective(self._handle, ctypes.byref(objective)))
        _check(self._lib, self._lib.inferno_get_iterations(self._handle, ctypes.byref(iters)))
        _check(self._lib, self._lib.inferno_get_checker_passed(self._handle, ctypes.byref(passed)))

        x, y = [], []
        if status.value == SolveStatus.OPTIMAL:
            nc, nr = self.num_cols, self.num_rows
            xbuf = (c_double * nc)()
            _check(self._lib, self._lib.inferno_get_solution(self._handle, xbuf, nc))
            x = list(xbuf)
            ybuf = (c_double * nr)()
            _check(self._lib, self._lib.inferno_get_duals(self._handle, ybuf, nr))
            y = list(ybuf)

        return Result(status.value, objective.value, iters.value, bool(passed.value), x, y)

    def close(self):
        if getattr(self, "_handle", None):
            self._lib.inferno_problem_destroy(self._handle)
            self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass  # interpreter teardown; nothing useful to do
