"""The contract between the language model and the solver.

Everything an LLM produces is untrusted input. This module is the only
place that decides whether a candidate model is allowed to reach Inferno,
and it decides deterministically: no network, no model, no heuristics that
could be talked out of a verdict. If `validate_model` returns errors, the
solver is not run. That rule has no exceptions and no override flag,
because the failure mode it prevents -- a plausible-looking but wrong
model solved to "optimality" and presented with real solver statistics --
is exactly the failure a judge or an operator cannot detect by eye.

The schema is deliberately narrow. It describes numbers, names and bounds
and nothing else: there is no expression language, no formula string, no
code field. An LLM cannot ask this schema to execute anything, because
there is nothing in it that executes.
"""

import json
import math
import re

# Mirrors the JSON Schema handed to the model for structured output. Kept
# next to the validator on purpose: when the schema and the validator
# disagree, the validator wins, and having both in one file makes the
# disagreement visible rather than silent.
MAX_VARIABLES = 2000
MAX_CONSTRAINTS = 2000
MAX_TERMS_PER_CONSTRAINT = 2000
MAX_ABS_COEFFICIENT = 1e12
MAX_ABS_BOUND = 1e15

NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_.\-]{0,63}$")


def _was_negative_infinity(value):
    """An explicit -inf lower bound means 'free below', which is different
    from an omitted lower bound, whose conventional default is 0."""
    return isinstance(value, (int, float)) and not isinstance(value, bool) \
        and math.isinf(float(value)) and float(value) < 0

VARIABLE_TYPES = ("continuous", "integer", "binary")
SENSES = ("minimize", "maximize")

JSON_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "required": ["name", "objective_sense", "variables", "objective", "constraints"],
    "properties": {
        "name": {"type": "string"},
        "objective_sense": {"type": "string", "enum": list(SENSES)},
        "objective_offset": {"type": "number"},
        "variables": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "required": ["name", "type", "lower_bound", "upper_bound"],
                "properties": {
                    "name": {"type": "string"},
                    "type": {"type": "string", "enum": list(VARIABLE_TYPES)},
                    "lower_bound": {"type": ["number", "null"]},
                    "upper_bound": {"type": ["number", "null"]},
                    "description": {"type": "string"},
                },
            },
        },
        "objective": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "required": ["variable", "coefficient"],
                "properties": {
                    "variable": {"type": "string"},
                    "coefficient": {"type": "number"},
                },
            },
        },
        "constraints": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "required": ["name", "terms", "lower_bound", "upper_bound"],
                "properties": {
                    "name": {"type": "string"},
                    "description": {"type": "string"},
                    "lower_bound": {"type": ["number", "null"]},
                    "upper_bound": {"type": ["number", "null"]},
                    "terms": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "additionalProperties": False,
                            "required": ["variable", "coefficient"],
                            "properties": {
                                "variable": {"type": "string"},
                                "coefficient": {"type": "number"},
                            },
                        },
                    },
                },
            },
        },
    },
}


class ValidationError(Exception):
    """Raised only by `require_valid`. Carries the full error list so the
    UI can show every problem at once instead of one per round trip."""

    def __init__(self, errors):
        self.errors = list(errors)
        super().__init__("; ".join(self.errors))


def _finite_number(value, where, errors, allow_none=False, limit=MAX_ABS_BOUND,
                   allow_infinite=False):
    """Numbers are where a malformed model does the most damage, because a
    NaN bound propagates silently through floating-point comparisons and a
    1e400 coefficient wrecks the factorization. Both are rejected here.

    `allow_infinite` is set for BOUND fields only, and it exists because
    validation has to be idempotent: this function normalizes an omitted
    bound to +/-inf, so a model that has already been through here must be
    accepted if it comes back -- which is exactly what happens every time
    the UI re-runs an edited model, and what a caller does when it passes a
    validated model straight to route_solve. Without this, validate(
    validate(m)) rejected its own output. Coefficients still refuse
    infinities outright: there is no such thing as an unbounded
    coefficient."""
    if value is None:
        if allow_none:
            return None
        errors.append(f"{where}: missing (null not allowed here)")
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        errors.append(f"{where}: must be a number, got {type(value).__name__}")
        return None
    v = float(value)
    if math.isnan(v):
        errors.append(f"{where}: is NaN")
        return None
    if math.isinf(v):
        if allow_infinite:
            return None  # same meaning as an omitted bound: unbounded on that side
        errors.append(f"{where}: is infinite (use null for unbounded)")
        return None
    if abs(v) > limit:
        errors.append(f"{where}: magnitude {v:g} exceeds the allowed limit {limit:g}")
        return None
    return v


def validate_model(raw):
    """Checks a candidate model and returns (normalized_model, errors).

    On any error the normalized model is None -- there is no partially
    accepted model, because a model that is 90% valid is still a model
    whose answer means nothing."""
    errors = []

    if isinstance(raw, (str, bytes)):
        try:
            raw = json.loads(raw)
        except (ValueError, TypeError) as exc:
            return None, [f"not valid JSON: {exc}"]

    if not isinstance(raw, dict):
        return None, [f"model must be a JSON object, got {type(raw).__name__}"]

    for field in ("objective_sense", "variables", "objective", "constraints"):
        if field not in raw:
            errors.append(f"missing required field '{field}'")
    if errors:
        return None, errors

    sense = raw.get("objective_sense")
    if sense not in SENSES:
        errors.append(f"objective_sense must be one of {SENSES}, got {sense!r}")

    name = raw.get("name") or "ai_model"
    if not isinstance(name, str):
        errors.append("name must be a string")
        name = "ai_model"

    offset = _finite_number(raw.get("objective_offset", 0.0), "objective_offset", errors,
                            allow_none=True)
    if offset is None:
        offset = 0.0

    # ---- variables ----
    variables = raw.get("variables")
    if not isinstance(variables, list):
        return None, errors + ["variables must be an array"]
    if not variables:
        errors.append("variables is empty: a model needs at least one variable")
    if len(variables) > MAX_VARIABLES:
        errors.append(f"too many variables ({len(variables)} > {MAX_VARIABLES})")

    index = {}
    norm_vars = []
    for i, v in enumerate(variables):
        where = f"variables[{i}]"
        if not isinstance(v, dict):
            errors.append(f"{where}: must be an object")
            continue
        vname = v.get("name")
        if not isinstance(vname, str) or not NAME_RE.match(vname):
            errors.append(
                f"{where}: name {vname!r} is not a valid identifier "
                "(letters, digits, _ . - ; must start with a letter or _)")
            continue
        if vname in index:
            errors.append(f"{where}: duplicate variable name {vname!r}")
            continue
        vtype = v.get("type")
        if vtype not in VARIABLE_TYPES:
            errors.append(f"{where} ({vname}): type must be one of {VARIABLE_TYPES}, got {vtype!r}")
            continue

        lo = _finite_number(v.get("lower_bound"), f"{where} ({vname}).lower_bound", errors,
                            allow_none=True, allow_infinite=True)
        hi = _finite_number(v.get("upper_bound"), f"{where} ({vname}).upper_bound", errors,
                            allow_none=True, allow_infinite=True)

        if vtype == "binary":
            # A binary variable whose box is not {0,1} is not binary. Rather
            # than silently "fixing" it -- which would change the model the
            # user was shown -- anything outside [0,1] is an error, and the
            # common case of omitted bounds is filled in.
            if lo is None:
                lo = 0.0
            if hi is None:
                hi = 1.0
            if lo < 0.0 or hi > 1.0:
                errors.append(
                    f"{where} ({vname}): binary variable must have bounds inside [0,1], "
                    f"got [{lo:g}, {hi:g}]")
                continue
        else:
            if lo is None:
                lo = -math.inf if _was_negative_infinity(v.get("lower_bound")) else 0.0
            if hi is None:
                hi = math.inf

        if lo > hi:
            errors.append(f"{where} ({vname}): lower_bound {lo:g} exceeds upper_bound {hi:g}")
            continue
        if vtype in ("integer", "binary"):
            # An integer column whose box contains no integer is infeasible
            # by construction, and the solver would spend a search proving
            # it. Say so up front instead.
            if math.floor(hi) < math.ceil(lo):
                errors.append(
                    f"{where} ({vname}): no integer lies in [{lo:g}, {hi:g}]")
                continue

        desc = v.get("description")
        if desc is not None and not isinstance(desc, str):
            errors.append(f"{where} ({vname}): description must be a string")
            desc = None

        index[vname] = len(norm_vars)
        norm_vars.append({
            "name": vname,
            "type": vtype,
            "lower_bound": lo,
            "upper_bound": hi,
            "description": desc or "",
        })

    # ---- objective ----
    objective = raw.get("objective")
    if not isinstance(objective, list):
        errors.append("objective must be an array of {variable, coefficient}")
        objective = []
    obj_coeff = {}
    for i, t in enumerate(objective):
        where = f"objective[{i}]"
        if not isinstance(t, dict):
            errors.append(f"{where}: must be an object")
            continue
        vname = t.get("variable")
        if vname not in index:
            errors.append(f"{where}: refers to unknown variable {vname!r}")
            continue
        if vname in obj_coeff:
            errors.append(f"{where}: duplicate objective term for {vname!r}")
            continue
        c = _finite_number(t.get("coefficient"), f"{where}.coefficient", errors,
                           limit=MAX_ABS_COEFFICIENT)
        if c is None:
            continue
        obj_coeff[vname] = c

    if not obj_coeff and not errors:
        errors.append("objective has no usable terms: nothing to optimize")

    # ---- constraints ----
    constraints = raw.get("constraints")
    if not isinstance(constraints, list):
        errors.append("constraints must be an array")
        constraints = []
    if len(constraints) > MAX_CONSTRAINTS:
        errors.append(f"too many constraints ({len(constraints)} > {MAX_CONSTRAINTS})")

    norm_cons = []
    seen_con_names = set()
    for i, c in enumerate(constraints):
        where = f"constraints[{i}]"
        if not isinstance(c, dict):
            errors.append(f"{where}: must be an object")
            continue
        cname = c.get("name")
        if not isinstance(cname, str) or not NAME_RE.match(cname):
            errors.append(f"{where}: name {cname!r} is not a valid identifier")
            continue
        if cname in seen_con_names:
            errors.append(f"{where}: duplicate constraint name {cname!r}")
            continue

        terms = c.get("terms")
        if not isinstance(terms, list):
            errors.append(f"{where} ({cname}): terms must be an array")
            continue
        if not terms:
            errors.append(f"{where} ({cname}): has no terms")
            continue
        if len(terms) > MAX_TERMS_PER_CONSTRAINT:
            errors.append(f"{where} ({cname}): too many terms ({len(terms)})")
            continue

        coeffs = {}
        bad_terms = False
        for k, t in enumerate(terms):
            tw = f"{where} ({cname}).terms[{k}]"
            if not isinstance(t, dict):
                errors.append(f"{tw}: must be an object")
                bad_terms = True
                continue
            vname = t.get("variable")
            if vname not in index:
                errors.append(f"{tw}: refers to unknown variable {vname!r}")
                bad_terms = True
                continue
            if vname in coeffs:
                errors.append(f"{tw}: duplicate term for {vname!r} in the same constraint")
                bad_terms = True
                continue
            coeff = _finite_number(t.get("coefficient"), f"{tw}.coefficient", errors,
                                   limit=MAX_ABS_COEFFICIENT)
            if coeff is None:
                bad_terms = True
                continue
            coeffs[vname] = coeff
        if bad_terms:
            continue

        lo = _finite_number(c.get("lower_bound"), f"{where} ({cname}).lower_bound", errors,
                            allow_none=True, allow_infinite=True)
        hi = _finite_number(c.get("upper_bound"), f"{where} ({cname}).upper_bound", errors,
                            allow_none=True, allow_infinite=True)
        if lo is None and hi is None:
            errors.append(
                f"{where} ({cname}): both bounds are null, so the row constrains nothing")
            continue
        lo = -math.inf if lo is None else lo
        hi = math.inf if hi is None else hi
        if lo > hi:
            errors.append(
                f"{where} ({cname}): lower_bound {lo:g} exceeds upper_bound {hi:g}")
            continue

        cdesc = c.get("description")
        if cdesc is not None and not isinstance(cdesc, str):
            cdesc = None

        seen_con_names.add(cname)
        norm_cons.append({
            "name": cname,
            "description": cdesc or "",
            "terms": [{"variable": v, "coefficient": coeffs[v]} for v in coeffs],
            "lower_bound": lo,
            "upper_bound": hi,
        })

    if errors:
        return None, errors

    model = {
        "name": name,
        "objective_sense": sense,
        "objective_offset": offset,
        "variables": norm_vars,
        "objective": [{"variable": v, "coefficient": obj_coeff[v]} for v in obj_coeff],
        "constraints": norm_cons,
    }
    return model, []


def require_valid(raw):
    model, errors = validate_model(raw)
    if errors:
        raise ValidationError(errors)
    return model


def format_model(model):
    """A human-readable rendering of a validated model. Shown to the user
    BEFORE the solver runs, so the formulation can be checked by a person
    rather than taken on faith from the language model that wrote it."""
    lines = []
    sense = "Minimize" if model["objective_sense"] == "minimize" else "Maximize"
    terms = " + ".join(
        f"{t['coefficient']:g}·{t['variable']}" for t in model["objective"])
    terms = terms.replace("+ -", "- ")
    lines.append(f"{sense}   {terms}" + (
        f" + {model['objective_offset']:g}" if model.get("objective_offset") else ""))
    lines.append("")
    lines.append("subject to")
    for c in model["constraints"]:
        expr = " + ".join(f"{t['coefficient']:g}·{t['variable']}" for t in c["terms"])
        expr = expr.replace("+ -", "- ")
        lo, hi = c["lower_bound"], c["upper_bound"]
        if lo == hi:
            body = f"{expr} = {lo:g}"
        elif math.isinf(lo):
            body = f"{expr} <= {hi:g}"
        elif math.isinf(hi):
            body = f"{expr} >= {lo:g}"
        else:
            body = f"{lo:g} <= {expr} <= {hi:g}"
        lines.append(f"  {c['name']}:  {body}")
    lines.append("")
    lines.append("bounds")
    for v in model["variables"]:
        if v["type"] == "binary":
            lines.append(f"  {v['name']} in {{0, 1}}")
        else:
            hi = "inf" if math.isinf(v["upper_bound"]) else f"{v['upper_bound']:g}"
            kind = " (integer)" if v["type"] == "integer" else ""
            lines.append(f"  {v['lower_bound']:g} <= {v['name']} <= {hi}{kind}")
    return "\n".join(lines)
