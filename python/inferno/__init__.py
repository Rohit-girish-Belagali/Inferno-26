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
from ctypes import POINTER, Structure, c_char_p, c_double, c_int, c_void_p

__all__ = ["Problem", "Result", "MipResult", "NodeEvent", "SolveStatus", "InfernoError",
           "load_library", "version"]

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


class NodeOutcome:
    ROOT = 0
    BRANCHED = 1
    INTEGER_FEASIBLE = 2
    INFEASIBLE = 3
    DOMINATED = 4
    GAP_CUT = 5
    RELAXATION_FAILED = 6

    _NAMES = {0: "root", 1: "branched", 2: "integer", 3: "infeasible",
              4: "dominated", 5: "gap-cut", 6: "unresolved"}

    @classmethod
    def name(cls, value):
        return cls._NAMES.get(value, f"unknown({value})")


class _CNodeEvent(Structure):
    """Mirrors inferno_node_event exactly. Field order and types must match
    the C struct or every read is garbage — ctypes cannot check this."""
    _fields_ = [
        ("node_index", c_int),
        ("depth", c_int),
        ("branch_var", c_int),
        ("branch_value", c_double),
        ("node_bound", c_double),
        ("incumbent", c_double),
        ("best_bound", c_double),
        ("outcome", c_int),
        ("elapsed_seconds", c_double),
        ("open_nodes", c_int),
    ]


class NodeEvent:
    """One node of a real branch-and-bound search, as the solver reported
    it. Nothing here is computed on the Python side."""

    __slots__ = ("node_index", "depth", "branch_var", "branch_value", "node_bound",
                 "incumbent", "best_bound", "outcome", "elapsed_seconds", "open_nodes")

    def __init__(self, c):
        self.node_index = c.node_index
        self.depth = c.depth
        self.branch_var = c.branch_var
        self.branch_value = c.branch_value
        self.node_bound = _from_sentinel(c.node_bound)
        self.incumbent = _from_sentinel(c.incumbent)
        self.best_bound = _from_sentinel(c.best_bound)
        self.outcome = c.outcome
        self.elapsed_seconds = c.elapsed_seconds
        self.open_nodes = c.open_nodes

    @property
    def outcome_name(self):
        return NodeOutcome.name(self.outcome)

    def as_dict(self):
        return {
            "node": self.node_index,
            "depth": self.depth,
            "branch_var": self.branch_var,
            "branch_value": self.branch_value,
            "node_bound": self.node_bound,
            "incumbent": self.incumbent,
            "best_bound": self.best_bound,
            "outcome": self.outcome_name,
            "elapsed": self.elapsed_seconds,
            "open_nodes": self.open_nodes,
        }


def _from_sentinel(v):
    """The C side represents 'no incumbent yet' / 'no bound yet' as the
    large finite INFINITY sentinel, because not every calling language can
    express IEEE infinity. Python can, so convert back rather than letting
    1e308 leak into a UI as a real number."""
    if v >= INFINITY:
        return float("inf")
    if v <= -INFINITY:
        return float("-inf")
    return v


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
    lib.inferno_problem_set_integer.restype = c_int
    lib.inferno_problem_set_integer.argtypes = [c_void_p, POINTER(c_int), c_int]
    lib.inferno_solve_mip.restype = c_int
    lib.inferno_solve_mip.argtypes = [c_void_p, c_double, c_int, c_double, c_int]
    for fn in ("inferno_get_nodes_explored", "inferno_get_proved_optimal",
               "inferno_get_event_count"):
        getattr(lib, fn).restype = c_int
        getattr(lib, fn).argtypes = [c_void_p, POINTER(c_int)]
    for fn in ("inferno_get_best_bound", "inferno_get_gap", "inferno_get_solve_seconds"):
        getattr(lib, fn).restype = c_int
        getattr(lib, fn).argtypes = [c_void_p, POINTER(c_double)]
    lib.inferno_get_events.restype = c_int
    lib.inferno_get_events.argtypes = [c_void_p, POINTER(_CNodeEvent), c_int, c_int,
                                       POINTER(c_int)]
    lib.inferno_node_outcome_string.restype = c_char_p
    lib.inferno_node_outcome_string.argtypes = [c_int]
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


class MipResult(Result):
    """A MILP answer. `proved_optimal` is the field that matters and is
    deliberately separate from `status`: branch-and-bound can return a
    perfectly good incumbent it could not prove optimal, and collapsing
    those two into one word is how a solver ends up overclaiming."""

    def __init__(self, status, objective, nodes, checker_passed, x, best_bound, gap,
                 proved_optimal, seconds):
        super().__init__(status, objective, nodes, checker_passed, x, [])
        self.nodes = nodes
        self.best_bound = best_bound
        self.gap = gap
        self.proved_optimal = proved_optimal
        self.seconds = seconds

    def __repr__(self):
        return (f"MipResult(status={self.status_name}, objective={self.objective!r}, "
                f"bound={self.best_bound!r}, gap={self.gap!r}, nodes={self.nodes}, "
                f"proved_optimal={self.proved_optimal}, checker_passed={self.checker_passed})")


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

    def set_integer(self, flags):
        """Marks which columns must be integral, making this a MILP. Pass
        None to clear every flag and go back to a pure LP."""
        if flags is None:
            _check(self._lib, self._lib.inferno_problem_set_integer(self._handle, None, 0))
            return
        buf = (c_int * len(flags))(*[1 if f else 0 for f in flags])
        _check(self._lib, self._lib.inferno_problem_set_integer(self._handle, buf, len(flags)))

    def solve_mip(self, time_limit=30.0, node_limit=200000, gap_tolerance=1e-6, presolve=True):
        """Runs branch-and-bound. This releases the GIL for the duration of
        the C call, so another thread may poll `events()` to watch the real
        search progress while this is running."""
        _check(self._lib, self._lib.inferno_solve_mip(
            self._handle, float(time_limit), int(node_limit), float(gap_tolerance),
            1 if presolve else 0))

        status, nodes, passed, proved = c_int(), c_int(), c_int(), c_int()
        objective, bound, gap, seconds = c_double(), c_double(), c_double(), c_double()
        _check(self._lib, self._lib.inferno_get_solve_status(self._handle, ctypes.byref(status)))
        _check(self._lib, self._lib.inferno_get_objective(self._handle, ctypes.byref(objective)))
        _check(self._lib, self._lib.inferno_get_nodes_explored(self._handle, ctypes.byref(nodes)))
        _check(self._lib, self._lib.inferno_get_checker_passed(self._handle, ctypes.byref(passed)))
        _check(self._lib, self._lib.inferno_get_best_bound(self._handle, ctypes.byref(bound)))
        _check(self._lib, self._lib.inferno_get_gap(self._handle, ctypes.byref(gap)))
        _check(self._lib, self._lib.inferno_get_proved_optimal(self._handle, ctypes.byref(proved)))
        _check(self._lib, self._lib.inferno_get_solve_seconds(self._handle, ctypes.byref(seconds)))

        x = []
        if status.value == SolveStatus.OPTIMAL:
            nc = self.num_cols
            xbuf = (c_double * nc)()
            _check(self._lib, self._lib.inferno_get_solution(self._handle, xbuf, nc))
            x = list(xbuf)

        return MipResult(status.value, _from_sentinel(objective.value), nodes.value,
                         bool(passed.value), x, _from_sentinel(bound.value),
                         _from_sentinel(gap.value), bool(proved.value), seconds.value)

    def event_count(self):
        """Node events recorded so far. Safe to call from another thread
        while solve_mip is running."""
        v = c_int()
        _check(self._lib, self._lib.inferno_get_event_count(self._handle, ctypes.byref(v)))
        return v.value

    def events(self, start=0, limit=1000):
        """Reads recorded node events. Safe to call from another thread
        while solve_mip is running — this is how live search progress
        reaches a UI without anything being simulated."""
        buf = (_CNodeEvent * limit)()
        written = c_int()
        _check(self._lib, self._lib.inferno_get_events(
            self._handle, buf, int(start), int(limit), ctypes.byref(written)))
        return [NodeEvent(buf[i]) for i in range(written.value)]

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
