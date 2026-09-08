"""OpenRouter client: the only place this application talks to a language
model, and the only place that needs an API key.

Three rules shape everything here.

  1. The key lives in the process environment and never leaves the server.
     No route echoes it, no log line prints it, and no response body
     contains it. The browser never sees it.
  2. Every failure is a value, not an exception that escapes. The caller
     gets an LlmResult with `ok=False` and a machine-readable `code`, and
     the application carries on. A language model being down is a degraded
     feature, never an outage.
  3. Retries are bounded and only for errors that could plausibly succeed
     on a second attempt. A 400 is not retried; a 502 is, twice.

The standard library does all the HTTP. That is deliberate: this project's
central claim is that its dependency graph is inspectable, and adding a
vendored HTTP stack for four requests would work against it.
"""

import json
import os
import ssl
import time
import urllib.error
import urllib.request

from model_schema import JSON_SCHEMA

OPENROUTER_URL = "https://openrouter.ai/api/v1/chat/completions"

DEFAULT_TEXT_MODEL = "google/gemini-2.5-flash"
DEFAULT_IMAGE_MODEL = "google/gemini-3.1-flash-image-preview"

# Error codes. These are the observability vocabulary from the brief; they
# appear in logs and in API responses so a failure during a demo can be
# named in one word instead of guessed at.
AI_DISABLED = "AI_DISABLED"
AI_TIMEOUT = "AI_TIMEOUT"
AI_RATE_LIMIT = "AI_RATE_LIMIT"
AI_HTTP_ERROR = "AI_HTTP_ERROR"
AI_NETWORK_ERROR = "AI_NETWORK_ERROR"
AI_INVALID_JSON = "AI_INVALID_JSON"
AI_EMPTY_RESPONSE = "AI_EMPTY_RESPONSE"
AI_INVALID_MODEL = "AI_INVALID_MODEL"
IMAGE_FAILURE = "IMAGE_FAILURE"

SYSTEM_PROMPT = """You are Inferno AI, an optimization copilot.

Your purpose is to help users formulate real-world optimization problems
and to explain solutions produced by the Inferno solver.

You are NOT the optimization solver. Inferno is a from-scratch LP/MILP/QP
engine and it is the only source of truth for any solver result.

Never invent:
- objective values
- bounds
- optimality gaps
- node counts
- runtimes
- feasibility or infeasibility verdicts
- optimality claims
- solver status

Your responsibilities:
1. Understand natural-language optimization problems.
2. Ask a clarifying question when a required number is genuinely missing.
3. Formulate valid LP/MILP models as structured JSON.
4. Explain what the variables and constraints mean.
5. Explain results that Inferno has actually produced.
6. Answer follow-up questions using the actual model and the actual result.

Rules that are absolute:
- Never return executable code, shell commands, SQL, or anything intended
  to be run. You return structured optimization models and prose only.
- Never claim a solution is optimal before Inferno has solved the model.
- When you are given Inferno's results, use those exact numbers. If a
  result contradicts what you expected, the result is right and you are
  wrong; say so plainly rather than adjusting the numbers.
- Linear models only. If a problem needs a product of two variables, say
  so and propose a linearization rather than writing a nonlinear term.

Modeling conventions:
- Every constraint is a row with a lower_bound and an upper_bound. Use
  null for an unbounded side. `x >= 5` is lower_bound 5, upper_bound null.
  An equality has both bounds equal.
- A fixed cost paid only when a facility runs needs a binary variable and
  a linking row: production - capacity * open <= 0.
- Give every variable and constraint a short description. It is what the
  user reads when checking your formulation."""

MODEL_INSTRUCTION = """Return a single JSON object describing the
optimization model. Use only these fields:

{
  "name": "short_identifier",
  "objective_sense": "minimize" | "maximize",
  "variables": [
    {"name": "ident", "type": "continuous"|"integer"|"binary",
     "lower_bound": number|null, "upper_bound": number|null,
     "description": "what this represents"}
  ],
  "objective": [{"variable": "ident", "coefficient": number}],
  "constraints": [
    {"name": "ident", "description": "what this enforces",
     "terms": [{"variable": "ident", "coefficient": number}],
     "lower_bound": number|null, "upper_bound": number|null}
  ]
}

Names must match /^[A-Za-z_][A-Za-z0-9_.-]{0,63}$/. Every variable
referenced in the objective or a constraint must be declared in
"variables". Do not repeat a variable within one constraint; combine the
coefficients instead. Return the JSON object and nothing else."""

CLARIFY_INSTRUCTION = """If the user's message is missing a number you
genuinely cannot proceed without, respond with a JSON object of the form
{"clarification": "your question"} instead of a model. Only do this when
the problem is truly underspecified -- prefer stating a clearly labelled
assumption and building the model."""


class LlmResult:
    """Success or failure, never a raised exception. `code` is one of the
    AI_* constants above when ok is False."""

    def __init__(self, ok, content=None, code=None, message=None, raw=None, attempts=1):
        self.ok = ok
        self.content = content
        self.code = code
        self.message = message
        self.raw = raw
        self.attempts = attempts

    def as_dict(self):
        return {"ok": self.ok, "code": self.code, "message": self.message,
                "attempts": self.attempts}


class Config:
    """Reads the three environment variables the brief specifies. Absence
    of a key is a supported, tested state -- not an error at startup."""

    def __init__(self, env=None):
        env = env if env is not None else os.environ
        self.api_key = (env.get("OPENROUTER_API_KEY") or "").strip()
        self.text_model = (env.get("OPENROUTER_TEXT_MODEL") or DEFAULT_TEXT_MODEL).strip()
        self.image_model = (env.get("OPENROUTER_IMAGE_MODEL") or DEFAULT_IMAGE_MODEL).strip()
        self.timeout = float(env.get("OPENROUTER_TIMEOUT_SECONDS") or 45.0)
        self.max_attempts = int(env.get("OPENROUTER_MAX_ATTEMPTS") or 3)
        self.referer = env.get("OPENROUTER_REFERER") or "https://github.com/team-inferno"
        self.title = env.get("OPENROUTER_TITLE") or "Inferno Optimization Copilot"

    @property
    def enabled(self):
        return bool(self.api_key)

    def public(self):
        """What the browser is allowed to know: whether AI is available and
        which models are configured. Never the key, not even a prefix."""
        return {
            "enabled": self.enabled,
            "text_model": self.text_model if self.enabled else None,
            "image_model": self.image_model if self.enabled else None,
            "reason": None if self.enabled else "OPENROUTER_API_KEY not configured",
        }


def _retryable(status):
    """Which HTTP statuses are worth a second attempt.

    5xx and 429 are the obvious ones. 404 is here for a reason found the
    hard way against the live gateway: OpenRouter answers 404 when it has
    no provider endpoint available for the requested model right now, not
    only when the model does not exist. The identical request succeeded on
    the next attempt seconds later. Treating 404 as permanent meant a
    transient upstream hiccup surfaced as a hard failure mid-demo, with the
    fallback summary shown when a retry would have produced the real
    answer. A genuinely wrong model id still fails -- it just fails after
    the retry budget rather than immediately.
    """
    return status == 429 or status == 404 or status >= 500


def _post_json(url, payload, headers, timeout):
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=body, headers=headers, method="POST")
    ctx = ssl.create_default_context()
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
        return resp.status, resp.read()


class OpenRouterClient:
    """A thin, defensive wrapper. `transport` exists so the whole failure
    matrix can be tested without a network: the tests inject a callable
    with the same signature as `_post_json` and drive every branch."""

    def __init__(self, config=None, transport=None, logger=None, sleep=time.sleep):
        self.config = config or Config()
        self._transport = transport or _post_json
        self._log = logger or (lambda *a, **k: None)
        self._sleep = sleep

    # -- low level -------------------------------------------------------
    def _chat(self, messages, response_format=None, max_tokens=2400, temperature=0.2):
        cfg = self.config
        if not cfg.enabled:
            self._log(AI_DISABLED, "no API key configured")
            return LlmResult(False, code=AI_DISABLED,
                             message="OPENROUTER_API_KEY not configured")

        payload = {
            "model": cfg.text_model,
            "messages": messages,
            "max_tokens": max_tokens,
            "temperature": temperature,
        }
        if response_format is not None:
            payload["response_format"] = response_format

        headers = {
            "Authorization": f"Bearer {cfg.api_key}",
            "Content-Type": "application/json",
            "HTTP-Referer": cfg.referer,
            "X-Title": cfg.title,
        }

        last = None
        for attempt in range(1, max(1, cfg.max_attempts) + 1):
            self._log("AI_REQUEST", f"model={cfg.text_model} attempt={attempt}")
            try:
                status, raw = self._transport(OPENROUTER_URL, payload, headers, cfg.timeout)
            except urllib.error.HTTPError as exc:
                code = AI_RATE_LIMIT if exc.code == 429 else AI_HTTP_ERROR
                detail = f"HTTP {exc.code}"
                last = LlmResult(False, code=code, message=detail, attempts=attempt)
                self._log(code, detail)
                if _retryable(exc.code):
                    if attempt < cfg.max_attempts:
                        self._sleep(min(2.0 * attempt, 4.0))
                        continue
                return last
            except (TimeoutError, urllib.error.URLError, OSError) as exc:
                is_timeout = isinstance(exc, TimeoutError) or "timed out" in str(exc).lower()
                code = AI_TIMEOUT if is_timeout else AI_NETWORK_ERROR
                last = LlmResult(False, code=code, message=str(exc), attempts=attempt)
                self._log(code, str(exc))
                if attempt < cfg.max_attempts:
                    self._sleep(min(1.0 * attempt, 3.0))
                    continue
                return last
            except Exception as exc:  # a transport that misbehaves is still just a failure
                last = LlmResult(False, code=AI_NETWORK_ERROR, message=str(exc), attempts=attempt)
                self._log(AI_NETWORK_ERROR, str(exc))
                return last

            if status == 429:
                self._log(AI_RATE_LIMIT, "HTTP 429")
                last = LlmResult(False, code=AI_RATE_LIMIT, message="rate limited",
                                 attempts=attempt)
                if attempt < cfg.max_attempts:
                    self._sleep(min(2.0 * attempt, 4.0))
                    continue
                return last
            if _retryable(status):
                self._log(AI_HTTP_ERROR, f"HTTP {status}")
                last = LlmResult(False, code=AI_HTTP_ERROR, message=f"HTTP {status}",
                                 attempts=attempt)
                if attempt < cfg.max_attempts:
                    self._sleep(min(2.0 * attempt, 4.0))
                    continue
                return last
            if status >= 400:
                self._log(AI_HTTP_ERROR, f"HTTP {status}")
                return LlmResult(False, code=AI_HTTP_ERROR, message=f"HTTP {status}",
                                 attempts=attempt)

            try:
                data = json.loads(raw.decode("utf-8") if isinstance(raw, bytes) else raw)
            except (ValueError, UnicodeDecodeError) as exc:
                self._log(AI_INVALID_JSON, f"response envelope: {exc}")
                return LlmResult(False, code=AI_INVALID_JSON,
                                 message="the AI gateway returned a malformed response",
                                 attempts=attempt)

            content = _extract_content(data)
            if not content:
                self._log(AI_EMPTY_RESPONSE, "no message content")
                return LlmResult(False, code=AI_EMPTY_RESPONSE,
                                 message="the AI returned an empty response", attempts=attempt)

            self._log("AI_SUCCESS", f"chars={len(content)}")
            return LlmResult(True, content=content, raw=data, attempts=attempt)

        return last or LlmResult(False, code=AI_NETWORK_ERROR, message="no attempt was made")

    # -- high level ------------------------------------------------------
    def formulate(self, history, user_message, current_model=None):
        """Natural language in, a candidate structured model out. The model
        is NOT validated here -- that is model_schema's job, and keeping the
        two separate is what guarantees the validator cannot be skipped by a
        client that talks to this class directly.

        `current_model` is what makes a what-if follow-up work. "What if
        demand rises 20%?" is meaningless without the model it refers to:
        with only the conversation summary in context the model would
        invent a fresh set of costs and the comparison against the previous
        run would be against a different problem entirely."""
        messages = [{"role": "system", "content": SYSTEM_PROMPT + "\n\n" + MODEL_INSTRUCTION
                     + "\n\n" + CLARIFY_INSTRUCTION}]
        if current_model:
            messages.append({
                "role": "system",
                "content": "The user already has this model loaded. If their message asks "
                           "for a change to it (a what-if, a different capacity, an added "
                           "constraint), return the COMPLETE modified model with every "
                           "unchanged variable, coefficient and bound carried over exactly "
                           "as they are here. Change only what was asked for. If they are "
                           "describing a genuinely new problem instead, ignore this "
                           "model.\n\n```json\n"
                           + json.dumps(current_model, default=_json_default)[:40000]
                           + "\n```"})
        messages.extend(_trim_history(history))
        messages.append({"role": "user", "content": user_message})

        response_format = {
            "type": "json_schema",
            "json_schema": {"name": "optimization_model", "strict": False,
                            "schema": _schema_with_clarification()},
        }
        res = self._chat(messages, response_format=response_format)
        if not res.ok:
            return res, None

        parsed = _parse_json_object(res.content)
        if parsed is None:
            self._log(AI_INVALID_JSON, "model payload was not a JSON object")
            return LlmResult(False, code=AI_INVALID_JSON,
                             message="the AI did not return valid JSON", raw=res.content), None
        return res, parsed

    def explain(self, model, result, question=None, history=None):
        """Explanation only. The result dict is passed verbatim and the
        prompt forbids changing any number in it."""
        payload = {
            "model": {
                "name": model.get("name"),
                "objective_sense": model.get("objective_sense"),
                "variables": model.get("variables"),
                "objective": model.get("objective"),
                "constraints": model.get("constraints"),
            },
            "inferno_result": result,
        }
        instruction = (
            "Below is the model Inferno solved and the result Inferno produced. "
            "Explain the result to a planner in plain language: what was chosen, why "
            "it makes sense given the costs and constraints, which constraints are "
            "binding, and what the trade-offs are.\n\n"
            "Use the exact numbers from inferno_result. Do not round them into "
            "different values, do not add statistics that are not there, and do not "
            "claim optimality unless proved_optimal is true. If proved_optimal is "
            "false, say the search found this solution but did not prove it optimal, "
            "and give the gap Inferno reported.\n\n"
            "Answer in markdown. Be concise."
        )
        if question:
            instruction += f"\n\nThe user's specific question: {question}"

        messages = [{"role": "system", "content": SYSTEM_PROMPT}]
        messages.extend(_trim_history(history or []))
        messages.append({"role": "user",
                         "content": instruction + "\n\n```json\n"
                         + json.dumps(payload, default=_json_default)[:60000] + "\n```"})
        return self._chat(messages, max_tokens=1200, temperature=0.3)

    def chat(self, history, user_message, context=None):
        """Free-form conversational turn — used for follow-ups that are not
        a new formulation, like 'which variables are binary?'."""
        messages = [{"role": "system", "content": SYSTEM_PROMPT}]
        if context:
            messages.append({
                "role": "system",
                "content": "Current session context (the model and the last Inferno result). "
                           "Use these exact numbers; never invent others.\n\n```json\n"
                           + json.dumps(context, default=_json_default)[:60000] + "\n```"})
        messages.extend(_trim_history(history))
        messages.append({"role": "user", "content": user_message})
        return self._chat(messages, max_tokens=1200, temperature=0.3)

    def image(self, prompt):
        """Optional illustration. Never on the critical path: every caller
        treats a failure here as 'no picture', not as an error."""
        cfg = self.config
        if not cfg.enabled:
            return LlmResult(False, code=AI_DISABLED, message="OPENROUTER_API_KEY not configured")
        payload = {
            "model": cfg.image_model,
            "messages": [{"role": "user", "content": prompt}],
            "modalities": ["image", "text"],
        }
        headers = {
            "Authorization": f"Bearer {cfg.api_key}",
            "Content-Type": "application/json",
            "HTTP-Referer": cfg.referer,
            "X-Title": cfg.title,
        }
        try:
            status, raw = self._transport(OPENROUTER_URL, payload, headers, cfg.timeout)
            if status >= 400:
                self._log(IMAGE_FAILURE, f"HTTP {status}")
                return LlmResult(False, code=IMAGE_FAILURE, message=f"HTTP {status}")
            data = json.loads(raw.decode("utf-8") if isinstance(raw, bytes) else raw)
        except Exception as exc:
            self._log(IMAGE_FAILURE, str(exc))
            return LlmResult(False, code=IMAGE_FAILURE, message=str(exc))

        url = _extract_image(data)
        if not url:
            self._log(IMAGE_FAILURE, "no image in response")
            return LlmResult(False, code=IMAGE_FAILURE, message="image unavailable")
        return LlmResult(True, content=url, raw=None)


# -- helpers -------------------------------------------------------------

def _json_default(o):
    return str(o)


def _schema_with_clarification():
    schema = json.loads(json.dumps(JSON_SCHEMA))
    # The model may legitimately answer with a question instead of a model.
    # Allowing that in the schema is what keeps a clarification from being
    # emitted as a malformed model.
    schema["additionalProperties"] = False
    schema["properties"]["clarification"] = {"type": "string"}
    schema["required"] = []
    return schema


def _trim_history(history, max_turns=12, max_chars=6000):
    """Bounded context. An unbounded transcript is both a cost problem and a
    reliability problem: the request eventually exceeds the model's window
    and every later turn fails."""
    out = []
    for turn in list(history)[-max_turns:]:
        role = turn.get("role")
        content = turn.get("content")
        if role not in ("user", "assistant") or not isinstance(content, str):
            continue
        out.append({"role": role, "content": content[:max_chars]})
    return out


def _extract_content(data):
    try:
        choices = data.get("choices") or []
        if not choices:
            return None
        message = choices[0].get("message") or {}
        content = message.get("content")
        if isinstance(content, list):  # some gateways return content parts
            parts = [p.get("text", "") for p in content if isinstance(p, dict)]
            content = "".join(parts)
        return content if isinstance(content, str) and content.strip() else None
    except (AttributeError, IndexError, TypeError):
        return None


def _extract_image(data):
    try:
        message = (data.get("choices") or [{}])[0].get("message") or {}
        for img in message.get("images") or []:
            url = (img.get("image_url") or {}).get("url")
            if url:
                return url
    except (AttributeError, IndexError, TypeError):
        return None
    return None


def _parse_json_object(text):
    """Models wrap JSON in prose and fences more often than they should.
    Peeling those off is fine; guessing at malformed JSON is not, so this
    only ever parses, never repairs."""
    if not isinstance(text, str):
        return None
    s = text.strip()
    if s.startswith("```"):
        s = s.split("\n", 1)[-1] if "\n" in s else s
        if s.endswith("```"):
            s = s[: s.rfind("```")]
        s = s.strip()
        if s.startswith("json"):
            s = s[4:].strip()
    try:
        obj = json.loads(s)
        return obj if isinstance(obj, dict) else None
    except ValueError:
        pass
    start, end = s.find("{"), s.rfind("}")
    if start >= 0 and end > start:
        try:
            obj = json.loads(s[start:end + 1])
            return obj if isinstance(obj, dict) else None
        except ValueError:
            return None
    return None
