"""The Inferno copilot server.

Serves the dashboard and a small JSON API. The whole design follows one
rule from the brief, and it is worth stating before the routes:

    Inferno is the source of truth. The language model is an interface.

Structurally that means the AI-touching routes (/api/formulate, /api/chat,
/api/explain, /api/image) can every one of them fail without affecting
/api/solve, /api/examples or /api/validate. There is no code path where a
solve waits on, depends on, or is influenced by a language model call, and
`test_failures.py` asserts exactly that.

Standard library only, for the same reason the rest of this project has no
third-party dependencies: the point is that you can read everything that
runs.
"""

import argparse
import collections
import json
import math
import os
import sys
import threading
import time
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import bridge  # noqa: E402
import examples  # noqa: E402
import llm  # noqa: E402
from model_schema import format_model, validate_model  # noqa: E402

STATIC_DIR = os.path.join(ROOT, "dashboard")
MAX_BODY = 4 * 1024 * 1024
MAX_RUNS = 64


class EventLog:
    """A bounded in-memory log with a hard rule: values that could contain a
    secret never get here. Every call site passes a code and a short
    message it constructed itself, so an API key cannot arrive by
    accident."""

    def __init__(self, capacity=400):
        self._entries = collections.deque(maxlen=capacity)
        self._lock = threading.Lock()

    def log(self, code, message=""):
        entry = {"t": time.time(), "code": str(code), "message": str(message)[:400]}
        with self._lock:
            self._entries.append(entry)
        print(f"[{time.strftime('%H:%M:%S')}] {entry['code']} {entry['message']}", flush=True)

    def recent(self, limit=100):
        with self._lock:
            return list(self._entries)[-limit:]


class Sessions:
    """Conversation state, keyed by a session id the browser generates.
    Held in memory only: this is a demo console, not a multi-tenant service,
    and persisting user problem descriptions to disk would be a privacy
    cost with no benefit here."""

    def __init__(self):
        self._data = {}
        self._lock = threading.Lock()

    def get(self, sid):
        with self._lock:
            return self._data.setdefault(sid, {"history": [], "model": None, "result": None,
                                               "previous": None})

    def update(self, sid, **kwargs):
        with self._lock:
            state = self._data.setdefault(sid, {"history": [], "model": None, "result": None,
                                                "previous": None})
            state.update(kwargs)
            return state


class App:
    def __init__(self, config=None, client=None):
        self.event_log = EventLog()
        # `log` is the bound method, not the object: every call site in this
        # class and in the client logs by calling app.log(code, message).
        self.log = self.event_log.log
        self.config = config or llm.Config()
        self.client = client or llm.OpenRouterClient(self.config, logger=self.log)
        self.sessions = Sessions()
        self.runs = collections.OrderedDict()
        self._runs_lock = threading.Lock()

    # -- run registry ----------------------------------------------------
    def start_run(self, model, time_limit, node_limit, gap_tolerance):
        run_id = uuid.uuid4().hex[:12]
        run = bridge.SolveRun(run_id, model, time_limit, node_limit, gap_tolerance)
        with self._runs_lock:
            self.runs[run_id] = run
            while len(self.runs) > MAX_RUNS:
                self.runs.popitem(last=False)
        run.start()
        return run

    def get_run(self, run_id):
        with self._runs_lock:
            return self.runs.get(run_id)

    # -- routes ----------------------------------------------------------
    def health(self):
        return {
            "ok": True,
            "solver": {"version": bridge.inferno.version(), "engine": "Inferno"},
            "ai": self.config.public(),
            "examples": len(examples.CATALOG),
        }

    def route_examples(self):
        return {"examples": examples.catalog()}

    def route_example_model(self, example_id):
        raw = examples.build(example_id)
        if raw is None:
            return {"error": f"no example named {example_id!r}"}, 404
        model, errors = validate_model(raw)
        if errors:
            # A built-in example failing its own validator is a bug in this
            # repository, not user error, and it should be loud.
            self.log("EXAMPLE_INVALID", f"{example_id}: {errors[0]}")
            return {"error": "built-in example failed validation", "details": errors}, 500
        return {"model": model, "preview": format_model(model), "source": "example"}

    def route_validate(self, body):
        model, errors = validate_model(body.get("model"))
        if errors:
            return {"valid": False, "errors": errors}
        return {"valid": True, "model": model, "preview": format_model(model)}

    def route_formulate(self, body):
        """Natural language -> candidate model -> deterministic validation.

        The order matters and is the whole safety story: an AI response that
        does not validate is returned to the UI as an error, and no solve is
        started. There is no branch here that reaches the solver with an
        unvalidated model."""
        message = str(body.get("message") or "").strip()
        sid = str(body.get("session") or "default")[:64]
        if not message:
            return {"error": "empty message"}, 400
        if len(message) > 8000:
            return {"error": "message too long"}, 400

        state = self.sessions.get(sid)
        result, parsed = self.client.formulate(state["history"], message,
                                               current_model=state.get("model"))
        if not result.ok:
            self.log("AI_FALLBACK", f"formulate failed: {result.code}")
            return {"ok": False, "ai": result.as_dict(),
                    "fallback": True,
                    "message": _friendly(result.code, result.message)}

        if isinstance(parsed, dict) and parsed.get("clarification") and not parsed.get("variables"):
            state["history"].append({"role": "user", "content": message})
            state["history"].append({"role": "assistant", "content": parsed["clarification"]})
            return {"ok": True, "mode": "clarification",
                    "clarification": parsed["clarification"], "ai": result.as_dict()}

        model, errors = validate_model(parsed)
        if errors:
            self.log(llm.AI_INVALID_MODEL, f"{len(errors)} problem(s): {errors[0]}")
            return {"ok": False, "mode": "invalid_model", "errors": errors,
                    "raw": parsed, "ai": result.as_dict(),
                    "message": "The AI generated a model that could not be validated. "
                               "The solver was not executed."}

        state["history"].append({"role": "user", "content": message})
        state["history"].append({
            "role": "assistant",
            "content": f"Formulated model '{model['name']}' with "
                       f"{len(model['variables'])} variables and "
                       f"{len(model['constraints'])} constraints."})
        self.sessions.update(sid, model=model)
        return {"ok": True, "mode": "model", "model": model,
                "preview": format_model(model), "ai": result.as_dict()}

    def route_solve(self, body):
        """Validates again, then solves. Re-validating a model that arrived
        from the UI is not redundant: the UI lets the user edit the
        formulation, and an edited model is exactly as untrusted as a
        generated one."""
        model, errors = validate_model(body.get("model"))
        if errors:
            return {"ok": False, "errors": errors,
                    "message": "The model could not be validated. The solver was not "
                               "executed."}, 400
        time_limit = _clamp(body.get("time_limit", 20.0), 0.1, 300.0)
        node_limit = int(_clamp(body.get("node_limit", 200000), 1, 5000000))
        gap = _clamp(body.get("gap_tolerance", 1e-6), 0.0, 1.0)
        run = self.start_run(model, time_limit, node_limit, gap)
        self.log("SOLVE_START",
                 f"run={run.id} vars={len(model['variables'])} rows={len(model['constraints'])}")
        return {"ok": True, "run": run.id}

    def route_run(self, run_id, since):
        run = self.get_run(run_id)
        if run is None:
            return {"error": "unknown run"}, 404
        snap = run.snapshot(since)
        if run.done and snap.get("result") is not None and since == 0:
            pass
        return snap

    def route_explain(self, body):
        model, errors = validate_model(body.get("model"))
        if errors:
            return {"ok": False, "message": "model failed validation", "errors": errors}, 400
        result = body.get("result")
        if not isinstance(result, dict):
            return {"ok": False, "message": "no solver result supplied"}, 400
        sid = str(body.get("session") or "default")[:64]
        state = self.sessions.get(sid)
        self.sessions.update(sid, model=model, previous=state.get("result"), result=result)

        res = self.client.explain(model, result, question=body.get("question"),
                                  history=state["history"])
        if not res.ok:
            self.log("AI_FALLBACK", f"explain failed: {res.code}")
            return {"ok": False, "ai": res.as_dict(),
                    "message": _friendly(res.code, res.message),
                    "deterministic": _deterministic_summary(model, result)}
        state["history"].append({"role": "assistant", "content": res.content[:4000]})
        return {"ok": True, "explanation": res.content, "ai": res.as_dict()}

    def route_chat(self, body):
        message = str(body.get("message") or "").strip()
        if not message:
            return {"error": "empty message"}, 400
        sid = str(body.get("session") or "default")[:64]
        state = self.sessions.get(sid)
        context = {"model": state.get("model"), "last_result": state.get("result"),
                   "previous_result": state.get("previous")}
        res = self.client.chat(state["history"], message, context=context)
        if not res.ok:
            self.log("AI_FALLBACK", f"chat failed: {res.code}")
            return {"ok": False, "ai": res.as_dict(),
                    "message": _friendly(res.code, res.message)}
        state["history"].append({"role": "user", "content": message})
        state["history"].append({"role": "assistant", "content": res.content[:4000]})
        return {"ok": True, "reply": res.content, "ai": res.as_dict()}

    def route_image(self, body):
        """Illustrations only, and never on the critical path. A failure
        here returns ok=False with a 200 so the UI treats it as 'no
        picture' rather than as a broken request."""
        prompt = str(body.get("prompt") or "").strip()
        if not prompt:
            return {"ok": False, "message": "image unavailable"}
        guard = ("A clean conceptual illustration for an industrial optimization console. "
                 "No text, no numbers, no charts, no graphs, no data labels. Subject: ")
        res = self.client.image(guard + prompt[:400])
        if not res.ok:
            return {"ok": False, "code": res.code, "message": "image unavailable"}
        return {"ok": True, "image": res.content}

    def route_logs(self, limit):
        return {"entries": self.event_log.recent(limit)}


def _clamp(v, lo, hi):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return lo
    if math.isnan(v):
        return lo
    return max(lo, min(hi, v))


def _friendly(code, message):
    return {
        llm.AI_DISABLED: "OPENROUTER_API_KEY not configured. Inferno solver is still available.",
        llm.AI_TIMEOUT: "AI request timed out. Inferno solver is still available.",
        llm.AI_RATE_LIMIT: "AI temporarily unavailable (rate limited). "
                           "Inferno solver is still available.",
        llm.AI_HTTP_ERROR: "AI gateway error. Inferno solver is still available.",
        llm.AI_NETWORK_ERROR: "AI unreachable. Inferno solver is still available.",
        llm.AI_INVALID_JSON: "The AI returned something that was not a valid model. "
                             "The solver was not executed.",
        llm.AI_EMPTY_RESPONSE: "The AI returned an empty response. "
                               "Inferno solver is still available.",
    }.get(code, message or "AI temporarily unavailable. Inferno solver is still available.")


def _deterministic_summary(model, result):
    """The explanation shown when the language model is unavailable.

    Written from the solver's own numbers with no narrative added, which is
    why it is safe to show unconditionally: it cannot be wrong about the
    result, only terse about it."""
    lines = []
    status = result.get("status")
    if status != "OPTIMAL":
        lines.append(f"Inferno returned status {status}.")
        return "\n".join(lines)
    sense = "Minimum" if model["objective_sense"] == "minimize" else "Maximum"
    lines.append(f"{sense} objective: {result.get('objective'):,.4f}")
    if result.get("is_mip"):
        bound = result.get("best_bound")
        label = "Best bound (upper)" if result.get("bound_is_upper") else "Best bound (lower)"
        if bound is not None:
            lines.append(f"{label}: {bound:,.4f}")
        gap = result.get("gap")
        if gap is not None:
            lines.append(f"Optimality gap: {gap * 100:.4f}%")
        lines.append(f"Nodes explored: {result.get('nodes')}")
        lines.append("Inferno proved this solution optimal."
                     if result.get("proved_optimal")
                     else "Inferno found this solution but did NOT prove it optimal.")
    else:
        lines.append(f"Simplex iterations: {result.get('iterations')}")
    lines.append("Independent checker: "
                 + ("PASS" if result.get("checker_passed") else "did not pass"))
    nonzero = [v for v in result.get("variables", []) if abs(v["value"]) > 1e-9]
    if nonzero:
        lines.append("")
        lines.append("Non-zero variables:")
        for v in nonzero[:40]:
            lines.append(f"  {v['name']} = {v['value']:,.6g}")
    binding = [r["name"] for r in result.get("rows", []) if r.get("binding")]
    if binding:
        lines.append("")
        lines.append("Binding constraints: " + ", ".join(binding[:20]))
    return "\n".join(lines)


CONTENT_TYPES = {".html": "text/html; charset=utf-8", ".css": "text/css; charset=utf-8",
                 ".js": "application/javascript; charset=utf-8", ".json": "application/json",
                 ".svg": "image/svg+xml", ".png": "image/png", ".csv": "text/csv"}


def make_handler(app):
    class Handler(BaseHTTPRequestHandler):
        server_version = "Inferno/0.6"
        protocol_version = "HTTP/1.1"

        def log_message(self, *args):
            pass  # the EventLog is the log; the default one is noise

        # -- helpers ------------------------------------------------------
        def _send(self, payload, status=200, ctype="application/json"):
            body = payload if isinstance(payload, bytes) else json.dumps(
                _sanitize(payload), allow_nan=False, default=_json_safe).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            try:
                self.wfile.write(body)
            except (BrokenPipeError, ConnectionResetError):
                pass  # the browser navigated away mid-response

        def _read_body(self):
            try:
                length = int(self.headers.get("Content-Length") or 0)
            except ValueError:
                return None
            if length <= 0 or length > MAX_BODY:
                return {} if length <= 0 else None
            raw = self.rfile.read(length)
            try:
                obj = json.loads(raw.decode("utf-8"))
            except (ValueError, UnicodeDecodeError):
                return None
            return obj if isinstance(obj, dict) else None

        def _dispatch(self, fn):
            """Every route runs inside this. An unhandled exception in a
            handler returns a 500 and keeps the server alive -- a demo that
            dies on one bad request is worse than one that reports the bad
            request."""
            try:
                out = fn()
            except Exception as exc:
                app.log("SERVER_ERROR", f"{type(exc).__name__}: {exc}")
                self._send({"error": "internal server error",
                            "detail": f"{type(exc).__name__}"}, 500)
                return
            if isinstance(out, tuple):
                payload, status = out
            else:
                payload, status = out, 200
            self._send(payload, status)

        # -- verbs --------------------------------------------------------
        def do_GET(self):
            parsed = urlparse(self.path)
            path, query = parsed.path, parse_qs(parsed.query)

            if path == "/api/health":
                return self._dispatch(app.health)
            if path == "/api/examples":
                return self._dispatch(app.route_examples)
            if path.startswith("/api/examples/"):
                eid = path[len("/api/examples/"):].strip("/")
                return self._dispatch(lambda: app.route_example_model(eid))
            if path.startswith("/api/solve/"):
                run_id = path[len("/api/solve/"):].strip("/")
                since = int((query.get("since") or ["0"])[0] or 0)
                return self._dispatch(lambda: app.route_run(run_id, max(0, since)))
            if path == "/api/logs":
                limit = int((query.get("limit") or ["100"])[0] or 100)
                return self._dispatch(lambda: app.route_logs(min(400, max(1, limit))))
            if path.startswith("/api/"):
                return self._send({"error": "unknown endpoint"}, 404)
            return self._serve_static(path)

        def do_POST(self):
            path = urlparse(self.path).path
            body = self._read_body()
            if body is None:
                return self._send({"error": "request body was not a JSON object"}, 400)
            routes = {
                "/api/formulate": app.route_formulate,
                "/api/validate": app.route_validate,
                "/api/solve": app.route_solve,
                "/api/explain": app.route_explain,
                "/api/chat": app.route_chat,
                "/api/image": app.route_image,
            }
            fn = routes.get(path)
            if fn is None:
                return self._send({"error": "unknown endpoint"}, 404)
            return self._dispatch(lambda: fn(body))

        # -- static -------------------------------------------------------
        def _serve_static(self, path):
            if path in ("/", "/index.html", "/copilot"):
                path = "/copilot.html"
            rel = path.lstrip("/")
            # Contain every request inside dashboard/. commonpath is the
            # check that matters: a crafted path like /../server/llm.py
            # must not be able to read source, let alone anything else.
            target = os.path.realpath(os.path.join(STATIC_DIR, rel))
            static_root = os.path.realpath(STATIC_DIR)
            if os.path.commonpath([target, static_root]) != static_root:
                return self._send({"error": "not found"}, 404)
            if not os.path.isfile(target):
                return self._send({"error": "not found"}, 404)
            ext = os.path.splitext(target)[1].lower()
            with open(target, "rb") as fh:
                data = fh.read()
            self._send(data, 200, CONTENT_TYPES.get(ext, "application/octet-stream"))

    return Handler


def _sanitize(obj):
    """Replaces every non-finite float with null, everywhere in the payload.

    This is not cosmetic. Python's json module emits bare `Infinity` for an
    IEEE infinity, which is not valid JSON: the browser's JSON.parse throws
    on it and the whole response is lost. It also broke the round trip --
    the validator normalizes an omitted bound to +inf internally, so a model
    that came back out of this API could not be sent back in. Converting to
    null at the boundary fixes both, and null is already what the schema
    means by 'unbounded'."""
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize(v) for v in obj]
    return obj


def _json_safe(o):
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return str(o)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Inferno optimization copilot server")
    parser.add_argument("--port", type=int, default=8420)
    parser.add_argument("--host", default="127.0.0.1")
    args = parser.parse_args(argv)

    app = App()
    httpd = ThreadingHTTPServer((args.host, args.port), make_handler(app))
    cfg = app.config.public()
    print("Inferno copilot server")
    print(f"  solver     : Inferno {bridge.inferno.version()}")
    print(f"  AI gateway : " + ("OpenRouter · " + cfg["text_model"] if cfg["enabled"]
                                else "DISABLED (OPENROUTER_API_KEY not configured)"))
    print(f"  offline    : {len(examples.CATALOG)} predefined models, solver fully usable")
    print(f"  listening  : http://{args.host}:{args.port}/")
    app.log("SERVER_START", f"ai_enabled={cfg['enabled']}")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
