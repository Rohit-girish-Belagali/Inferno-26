"""Validated model -> Inferno -> real results.

This is the only path from an AI-generated model into the solver, and it
starts after `model_schema.validate_model` has already accepted the model.
Nothing here re-interprets the user's intent; it compiles an already
checked structure into the sparse form the C ABI wants, runs the real
solver, and copies the real numbers back out.

Two things it deliberately does NOT do:

  - It never computes an objective, bound or gap itself. Every number the
    UI shows comes from `inferno_get_*`. The one arithmetic operation here
    is the sign flip for a maximization, which is documented below and is
    a property of the formulation, not an estimate.
  - It never falls back to an approximation when the solver fails. A
    failed solve is reported as a failed solve.
"""

import math
import os
import sys
import threading
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "python"))

import inferno  # noqa: E402

INFINITY = 1.0e308


def _to_api(v):
    """Inferno's C ABI represents unbounded as a large finite sentinel, so
    that callers in languages without IEEE infinity can still express it."""
    if v == math.inf:
        return INFINITY
    if v == -math.inf:
        return -INFINITY
    return float(v)


def compile_model(model):
    """Turns a validated model into CSC arrays plus the integrality flags.

    Column-major by construction because that is what the solver wants: it
    accumulates each variable's coefficients across all rows first, then
    lays them out contiguously, so no transpose is needed anywhere.
    """
    var_index = {v["name"]: i for i, v in enumerate(model["variables"])}
    n = len(model["variables"])
    m = len(model["constraints"])

    cols = [[] for _ in range(n)]
    for i, c in enumerate(model["constraints"]):
        for t in c["terms"]:
            cols[var_index[t["variable"]]].append((i, t["coefficient"]))

    col_ptr, row_idx, values = [0], [], []
    for j in range(n):
        for (i, coeff) in sorted(cols[j]):
            row_idx.append(i)
            values.append(coeff)
        col_ptr.append(len(row_idx))

    # Inferno minimizes. A maximization is solved as the minimization of
    # the negated objective, and the sign is put back when the result is
    # read. This is an exact restatement of the same problem, not a
    # transformation that loses anything -- but it does swap the meaning of
    # the bound, which `run_solve` handles explicitly rather than leaving
    # to be noticed later.
    sign = 1.0 if model["objective_sense"] == "minimize" else -1.0
    obj = [0.0] * n
    for t in model["objective"]:
        obj[var_index[t["variable"]]] = sign * t["coefficient"]

    col_lo = [_to_api(v["lower_bound"]) for v in model["variables"]]
    col_hi = [_to_api(v["upper_bound"]) for v in model["variables"]]
    row_lo = [_to_api(c["lower_bound"]) for c in model["constraints"]]
    row_hi = [_to_api(c["upper_bound"]) for c in model["constraints"]]
    is_integer = [1 if v["type"] in ("integer", "binary") else 0 for v in model["variables"]]

    return {
        "num_rows": m,
        "num_cols": n,
        "col_ptr": col_ptr,
        "row_idx": row_idx,
        "values": values,
        "obj": obj,
        "col_lo": col_lo,
        "col_hi": col_hi,
        "row_lo": row_lo,
        "row_hi": row_hi,
        "is_integer": is_integer,
        "sign": sign,
    }


class SolveRun:
    """One solve, running on its own thread so the HTTP server stays
    responsive and the UI can poll genuine progress.

    The events it serves are read out of the solver's own event log via the
    C ABI while the search is still running -- ctypes releases the GIL for
    the duration of the foreign call, so the polling thread really does
    observe a live search rather than a replay.
    """

    def __init__(self, run_id, model, time_limit=20.0, node_limit=200000, gap_tolerance=1e-6):
        self.id = run_id
        self.model = model
        self.time_limit = float(time_limit)
        self.node_limit = int(node_limit)
        self.gap_tolerance = float(gap_tolerance)
        self.state = "queued"
        self.error = None
        self.result = None
        self.started_at = None
        self.finished_at = None
        self._problem = None
        self._lock = threading.Lock()
        self._events = []
        self._thread = None
        self._compiled = compile_model(model)

    # -- lifecycle -------------------------------------------------------
    def start(self):
        self._thread = threading.Thread(target=self._run, name=f"solve-{self.id}", daemon=True)
        self._thread.start()
        return self

    def join(self, timeout=None):
        if self._thread is not None:
            self._thread.join(timeout)
        return self

    @property
    def done(self):
        return self.state in ("finished", "failed")

    def _run(self):
        c = self._compiled
        self.started_at = time.time()
        self.state = "solving"
        try:
            p = inferno.Problem.from_arrays(
                c["num_rows"], c["num_cols"], c["col_ptr"], c["row_idx"], c["values"],
                c["obj"], c["col_lo"], c["col_hi"], c["row_lo"], c["row_hi"])
            with self._lock:
                self._problem = p
            has_integers = any(c["is_integer"])
            if has_integers:
                p.set_integer(c["is_integer"])
                r = p.solve_mip(time_limit=self.time_limit, node_limit=self.node_limit,
                                gap_tolerance=self.gap_tolerance, presolve=True)
            else:
                r = p.solve(presolve=True)
            self.result = self._read_result(r, has_integers)
            self.state = "finished"
        except Exception as exc:  # the solver failing must not take the server with it
            self.error = f"{type(exc).__name__}: {exc}"
            self.state = "failed"
        finally:
            self.finished_at = time.time()
            # Drain the final events before the handle goes away.
            self.events()
            with self._lock:
                p = self._problem
                self._problem = None
            if p is not None:
                try:
                    p.close()
                except Exception:
                    pass

    # -- results ---------------------------------------------------------
    def _read_result(self, r, has_integers):
        c = self._compiled
        sign = c["sign"]
        names = [v["name"] for v in self.model["variables"]]
        offset = self.model.get("objective_offset", 0.0)

        # Sign restoration for a maximization. The solver minimized -f, so
        # the user-facing objective is -(what it returned). The bound flips
        # with it: a LOWER bound on the minimization of -f is an UPPER
        # bound on the maximization of f, which is what `bound_is_upper`
        # tells the UI so it can label the number correctly instead of
        # calling every bound a lower bound.
        def user_obj(v):
            if v is None or not math.isfinite(v):
                return None
            return sign * v + offset

        out = {
            "status": r.status_name,
            "checker_passed": bool(r.checker_passed),
            "objective": user_obj(r.objective),
            "objective_sense": self.model["objective_sense"],
            "bound_is_upper": sign < 0,
            "variables": [],
            "is_mip": has_integers,
            # Read at result time, not at thread-exit time: _read_result runs
            # before the `finally` block that stamps finished_at, so keying
            # off that field reported null for every solve.
            "seconds": round(time.time() - self.started_at, 4) if self.started_at else None,
        }
        if r.status_name == "OPTIMAL" and r.x:
            for i, name in enumerate(names):
                v = self.model["variables"][i]
                out["variables"].append({
                    "name": name,
                    "value": r.x[i],
                    "type": v["type"],
                    "description": v["description"],
                    "lower_bound": v["lower_bound"] if math.isfinite(v["lower_bound"]) else None,
                    "upper_bound": v["upper_bound"] if math.isfinite(v["upper_bound"]) else None,
                })
            out["rows"] = self._row_report(r.x)

        if has_integers:
            out.update({
                "best_bound": user_obj(r.best_bound),
                "gap": r.gap if math.isfinite(r.gap) else None,
                "nodes": r.nodes,
                "proved_optimal": bool(r.proved_optimal),
                "solver_seconds": round(r.seconds, 4),
            })
        else:
            out.update({
                "iterations": r.iterations,
                "duals": [
                    {"name": self.model["constraints"][i]["name"], "dual": r.y[i]}
                    for i in range(len(r.y))
                ] if r.y else [],
            })
        return out

    def _row_report(self, x):
        """Row activities and which constraints are binding, recomputed from
        the model the user was shown. This is a second, independent pass over
        the same data the solver used -- the point is that a binding-constraint
        claim in the UI can be traced to arithmetic, not to a narrative."""
        var_index = {v["name"]: i for i, v in enumerate(self.model["variables"])}
        rows = []
        for c in self.model["constraints"]:
            activity = sum(t["coefficient"] * x[var_index[t["variable"]]] for t in c["terms"])
            lo, hi = c["lower_bound"], c["upper_bound"]
            tol = 1e-7 * max(1.0, abs(activity))
            binding_lo = math.isfinite(lo) and abs(activity - lo) <= tol
            binding_hi = math.isfinite(hi) and abs(activity - hi) <= tol
            slack = None
            if math.isfinite(hi):
                slack = hi - activity
            elif math.isfinite(lo):
                slack = activity - lo
            rows.append({
                "name": c["name"],
                "description": c["description"],
                "activity": activity,
                "lower_bound": lo if math.isfinite(lo) else None,
                "upper_bound": hi if math.isfinite(hi) else None,
                "binding": bool(binding_lo or binding_hi),
                "slack": slack,
            })
        return rows

    # -- live events -----------------------------------------------------
    def events(self, start=0):
        """Returns every recorded node event from `start` onwards. Pulls any
        new ones out of the solver first, including while it is still
        running."""
        with self._lock:
            p = self._problem
        if p is not None:
            try:
                while True:
                    fetched = p.events(start=len(self._events), limit=512)
                    if not fetched:
                        break
                    self._events.extend(e.as_dict() for e in fetched)
                    if len(fetched) < 512:
                        break
            except Exception:
                # A read failing mid-solve costs live progress, nothing more.
                pass
        return self._events[start:]

    def snapshot(self, since=0):
        evs = self.events(since)
        return {
            "id": self.id,
            "state": self.state,
            "error": self.error,
            "since": since,
            "events": evs,
            "event_count": len(self._events),
            "result": self.result,
            "elapsed": round((self.finished_at or time.time()) - self.started_at, 3)
            if self.started_at else 0.0,
        }


def solve_sync(model, time_limit=20.0, node_limit=200000, gap_tolerance=1e-6):
    """Blocking convenience wrapper, used by tests and by the offline CLI."""
    run = SolveRun("sync", model, time_limit, node_limit, gap_tolerance).start()
    run.join()
    return run
