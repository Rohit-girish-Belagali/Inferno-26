"""The copilot test matrix.

Two halves, and the split is the point of the whole exercise.

The first half checks that a language model cannot get a bad model past
the validator: malformed JSON, unknown variable references, duplicate
names, NaN and infinite numbers, binary variables with non-binary bounds,
a missing objective. Every one of these must be REJECTED, and rejection
must mean the solver was never called.

The second half checks that the language model failing does not matter:
missing key, timeout, HTTP 500, rate limit, network error, empty
response, malformed content, image failure. After each one the solver must
still solve a real model and the independent checker must still pass.

Every AI call is driven through an injected transport, so the whole matrix
runs with no API key, no network and no cost -- which also means it runs in
CI on a machine that has neither.
"""

import json
import os
import sys
import unittest
import urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
SERVER = os.path.dirname(HERE)
sys.path.insert(0, SERVER)

import app as app_module  # noqa: E402
import bridge  # noqa: E402
import examples  # noqa: E402
import llm  # noqa: E402
from model_schema import validate_model  # noqa: E402


GOOD_MODEL = {
    "name": "plants",
    "objective_sense": "minimize",
    "variables": [
        {"name": "open_A", "type": "binary", "lower_bound": 0, "upper_bound": 1,
         "description": "Open plant A"},
        {"name": "prod_A", "type": "continuous", "lower_bound": 0, "upper_bound": 500,
         "description": "Units at A"},
    ],
    "objective": [{"variable": "open_A", "coefficient": 50000},
                  {"variable": "prod_A", "coefficient": 20}],
    "constraints": [
        {"name": "demand", "terms": [{"variable": "prod_A", "coefficient": 1}],
         "lower_bound": 300, "upper_bound": None, "description": "meet demand"},
        {"name": "link_A", "terms": [{"variable": "prod_A", "coefficient": 1},
                                     {"variable": "open_A", "coefficient": -500}],
         "lower_bound": None, "upper_bound": 0, "description": "link"},
    ],
}


def envelope(content):
    """An OpenRouter chat-completions response body."""
    return json.dumps({"choices": [{"message": {"content": content}}]}).encode()


def transport_returning(content, status=200):
    calls = []

    def t(url, payload, headers, timeout):
        calls.append({"url": url, "payload": payload, "headers": headers})
        return status, envelope(content) if isinstance(content, str) else content
    t.calls = calls
    return t


def transport_raising(exc):
    calls = []

    def t(url, payload, headers, timeout):
        calls.append(1)
        raise exc
    t.calls = calls
    return t


def make_client(transport, key="test-key", attempts=3):
    cfg = llm.Config({"OPENROUTER_API_KEY": key, "OPENROUTER_MAX_ATTEMPTS": str(attempts)})
    return llm.OpenRouterClient(cfg, transport=transport, sleep=lambda _s: None)


def make_app(client=None, key=""):
    cfg = llm.Config({"OPENROUTER_API_KEY": key})
    return app_module.App(config=cfg, client=client or llm.OpenRouterClient(
        cfg, transport=transport_raising(AssertionError("transport must not be reached"))))


def solve_now(application, model, time_limit=20.0):
    out = application.route_solve({"model": model, "time_limit": time_limit})
    assert out.get("ok"), out
    run = application.get_run(out["run"])
    run.join(60)
    return run


# ---------------------------------------------------------------- validation

class ValidatorMatrix(unittest.TestCase):
    """Invalid LLM output must NEVER reach the solver."""

    def assert_rejected(self, raw, needle):
        model, errors = validate_model(raw)
        self.assertIsNone(model, f"expected rejection, got a model for {needle}")
        self.assertTrue(errors)
        joined = " | ".join(errors).lower()
        self.assertIn(needle.lower(), joined, f"errors were: {errors}")

    def test_valid_model_accepted(self):
        model, errors = validate_model(GOOD_MODEL)
        self.assertEqual(errors, [])
        self.assertEqual(len(model["variables"]), 2)

    def test_malformed_json(self):
        self.assert_rejected("{not json at all", "not valid json")

    def test_not_an_object(self):
        self.assert_rejected("[1, 2, 3]", "must be a json object")

    def test_missing_objective_field(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        del bad["objective"]
        self.assert_rejected(bad, "missing required field 'objective'")

    def test_missing_variables_field(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        del bad["variables"]
        self.assert_rejected(bad, "missing required field 'variables'")

    def test_duplicate_variable(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["variables"].append(dict(bad["variables"][0]))
        self.assert_rejected(bad, "duplicate variable name")

    def test_unknown_variable_in_constraint(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["constraints"][0]["terms"][0]["variable"] = "ghost"
        self.assert_rejected(bad, "unknown variable")

    def test_unknown_variable_in_objective(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["objective"][0]["variable"] = "ghost"
        self.assert_rejected(bad, "unknown variable")

    def test_duplicate_term_in_one_constraint(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["constraints"][1]["terms"].append({"variable": "prod_A", "coefficient": 3})
        self.assert_rejected(bad, "duplicate term")

    def test_inverted_bounds(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["variables"][1]["lower_bound"] = 900
        self.assert_rejected(bad, "exceeds upper_bound")

    def test_binary_with_non_binary_bounds(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["variables"][0]["upper_bound"] = 7
        self.assert_rejected(bad, "inside [0,1]")

    def test_integer_box_with_no_integer(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["variables"][0]["type"] = "integer"
        bad["variables"][0]["lower_bound"] = 0.2
        bad["variables"][0]["upper_bound"] = 0.8
        self.assert_rejected(bad, "no integer lies in")

    def test_nan_coefficient(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["objective"][0]["coefficient"] = float("nan")
        self.assert_rejected(bad, "is nan")

    def test_infinite_coefficient(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["constraints"][0]["terms"][0]["coefficient"] = float("inf")
        self.assert_rejected(bad, "is infinite")

    def test_oversized_coefficient(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["objective"][0]["coefficient"] = 1e30
        self.assert_rejected(bad, "exceeds the allowed limit")

    def test_bad_variable_type(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["variables"][0]["type"] = "quantum"
        self.assert_rejected(bad, "type must be one of")

    def test_bad_objective_sense(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["objective_sense"] = "optimise a bit"
        self.assert_rejected(bad, "objective_sense must be one of")

    def test_illegal_variable_name(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["variables"][0]["name"] = "rm -rf /"
        self.assert_rejected(bad, "not a valid identifier")

    def test_constraint_with_no_bounds(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["constraints"][0]["lower_bound"] = None
        bad["constraints"][0]["upper_bound"] = None
        self.assert_rejected(bad, "constrains nothing")

    def test_empty_constraint_terms(self):
        bad = json.loads(json.dumps(GOOD_MODEL))
        bad["constraints"][0]["terms"] = []
        self.assert_rejected(bad, "has no terms")


# ------------------------------------------------------- AI failure handling

class AiFailureMatrix(unittest.TestCase):
    """Every one of these is a failure the demo could actually hit. After
    each, the solver must still work."""

    def assert_solver_still_works(self, application):
        model, errors = validate_model(examples.build("plant_selection"))
        self.assertEqual(errors, [])
        run = solve_now(application, model)
        self.assertEqual(run.state, "finished", run.error)
        self.assertEqual(run.result["status"], "OPTIMAL")
        self.assertTrue(run.result["checker_passed"])
        self.assertTrue(run.result["proved_optimal"])

    def test_missing_api_key(self):
        a = make_app(key="")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertFalse(out["ok"])
        self.assertEqual(out["ai"]["code"], llm.AI_DISABLED)
        self.assertIn("still available", out["message"])
        self.assert_solver_still_works(a)

    def test_timeout(self):
        t = transport_raising(TimeoutError("timed out"))
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertEqual(out["ai"]["code"], llm.AI_TIMEOUT)
        self.assertEqual(len(t.calls), 3, "a timeout should be retried, bounded at 3")
        self.assert_solver_still_works(a)

    def test_http_500_retries_then_gives_up(self):
        t = transport_returning("", status=500)
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertEqual(out["ai"]["code"], llm.AI_HTTP_ERROR)
        self.assertEqual(len(t.calls), 3)
        self.assert_solver_still_works(a)

    def test_rate_limit(self):
        t = transport_returning("", status=429)
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertEqual(out["ai"]["code"], llm.AI_RATE_LIMIT)
        self.assertIn("still available", out["message"])
        self.assert_solver_still_works(a)

    def test_http_404_is_retried(self):
        """OpenRouter answers 404 when it has no provider endpoint free for
        the model right now, not only when the model does not exist. Seen
        live: the same request failed once and succeeded seconds later, so
        treating 404 as permanent turned a transient hiccup into a visible
        mid-demo failure."""
        t = transport_returning("", status=404)
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertEqual(out["ai"]["code"], llm.AI_HTTP_ERROR)
        self.assertEqual(len(t.calls), 3)
        self.assert_solver_still_works(a)

    def test_404_that_recovers_on_retry_returns_the_model(self):
        calls = {"n": 0}

        def t(url, payload, headers, timeout):
            calls["n"] += 1
            if calls["n"] == 1:
                return 404, b"{}"
            return 200, envelope(json.dumps(GOOD_MODEL))

        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertTrue(out["ok"], out)
        self.assertEqual(out["mode"], "model")
        self.assertEqual(calls["n"], 2)

    def test_http_400_is_not_retried(self):
        t = transport_returning("", status=400)
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertEqual(out["ai"]["code"], llm.AI_HTTP_ERROR)
        self.assertEqual(len(t.calls), 1, "a 400 cannot become valid on a retry")

    def test_network_failure(self):
        t = transport_raising(urllib.error.URLError("name resolution failed"))
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertEqual(out["ai"]["code"], llm.AI_NETWORK_ERROR)
        self.assert_solver_still_works(a)

    def test_empty_response(self):
        t = transport_returning(json.dumps({"choices": []}).encode())
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertEqual(out["ai"]["code"], llm.AI_EMPTY_RESPONSE)
        self.assert_solver_still_works(a)

    def test_malformed_json_content(self):
        t = transport_returning("here is your model: {broken")
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertEqual(out["ai"]["code"], llm.AI_INVALID_JSON)
        self.assert_solver_still_works(a)

    def test_malformed_gateway_envelope(self):
        t = transport_returning(b"<html>502 Bad Gateway</html>")
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertEqual(out["ai"]["code"], llm.AI_INVALID_JSON)

    def test_invalid_model_never_reaches_the_solver(self):
        """The headline safety property, asserted directly: the AI returns
        well-formed JSON describing a broken model, and no run is created."""
        broken = json.loads(json.dumps(GOOD_MODEL))
        broken["constraints"][0]["terms"][0]["variable"] = "does_not_exist"
        t = transport_returning(json.dumps(broken))
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertFalse(out["ok"])
        self.assertEqual(out["mode"], "invalid_model")
        self.assertIn("solver was not executed", out["message"])
        self.assertEqual(len(a.runs), 0, "a rejected model must not start a solve")
        self.assert_solver_still_works(a)

    def test_image_failure_is_isolated(self):
        t = transport_returning("", status=500)
        a = make_app(client=make_client(t), key="k")
        out = a.route_image({"prompt": "a refinery at dusk"})
        self.assertFalse(out["ok"])
        self.assertEqual(out["message"], "image unavailable")
        self.assert_solver_still_works(a)

    def test_api_key_never_leaves_the_server(self):
        t = transport_returning(json.dumps(GOOD_MODEL))
        a = make_app(client=make_client(t, key="super-secret-key"), key="super-secret-key")
        out = a.route_formulate({"message": "pick a plant"})
        blob = json.dumps(out) + json.dumps(a.health()) + json.dumps(a.route_logs(200))
        self.assertNotIn("super-secret-key", blob)
        # and it IS actually being sent upstream, so the check above is meaningful
        self.assertIn("super-secret-key", t.calls[0]["headers"]["Authorization"])


# ------------------------------------------------------------ happy path

class SuccessPath(unittest.TestCase):

    def test_valid_ai_response_formulates_and_solves(self):
        t = transport_returning(json.dumps(examples.build("plant_selection")))
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "three plants, fixed costs, meet 1000 units"})
        self.assertTrue(out["ok"])
        self.assertEqual(out["mode"], "model")
        self.assertIn("Minimize", out["preview"])

        run = solve_now(a, out["model"])
        self.assertEqual(run.result["status"], "OPTIMAL")
        self.assertTrue(run.result["checker_passed"])
        self.assertTrue(run.result["proved_optimal"])
        self.assertAlmostEqual(run.result["objective"], 127200.0, places=6)

    def test_ai_fenced_json_is_unwrapped(self):
        fenced = "```json\n" + json.dumps(GOOD_MODEL) + "\n```"
        t = transport_returning(fenced)
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "one plant"})
        self.assertTrue(out["ok"], out)

    def test_clarification_instead_of_a_model(self):
        t = transport_returning(json.dumps({"clarification": "What is the demand?"}))
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "pick a plant"})
        self.assertTrue(out["ok"])
        self.assertEqual(out["mode"], "clarification")
        self.assertEqual(len(a.runs), 0)

    def test_explanation_receives_the_real_numbers(self):
        captured = {}

        def t(url, payload, headers, timeout):
            captured["payload"] = payload
            return 200, envelope("Inferno proved this optimal at 127,200.")
        a = make_app(client=make_client(t), key="k")
        model, _ = validate_model(examples.build("plant_selection"))
        run = solve_now(a, model)
        out = a.route_explain({"model": model, "result": run.result})
        self.assertTrue(out["ok"])
        sent = json.dumps(captured["payload"])
        self.assertIn("127200", sent.replace(",", ""))
        self.assertIn("proved_optimal", sent)
        self.assertIn("never invent", json.dumps(captured["payload"]).lower())

    def test_follow_up_scenario_compares_two_real_runs(self):
        """The what-if flow: solve, raise demand 20%, solve again. Both
        numbers must come from separate real solves."""
        a = make_app(key="")
        base, errors = validate_model(examples.build("plant_selection"))
        self.assertEqual(errors, [])
        first = solve_now(a, base)

        raised = json.loads(json.dumps(base))
        for c in raised["constraints"]:
            if c["name"] == "demand":
                c["lower_bound"] = c["lower_bound"] * 1.2
        raised, errors = validate_model(raised)
        self.assertEqual(errors, [])
        second = solve_now(a, raised)

        self.assertEqual(first.result["status"], "OPTIMAL")
        self.assertEqual(second.result["status"], "OPTIMAL")
        self.assertTrue(second.result["checker_passed"])
        self.assertGreater(second.result["objective"], first.result["objective"],
                           "more demand cannot be cheaper")

    def test_what_if_follow_up_carries_the_model_forward(self):
        """A what-if is meaningless without the model it refers to. The
        current model must be sent as context, or the AI re-derives a
        different problem and the comparison compares nothing."""
        captured = {}

        # Stateful on purpose: the first call is the original problem, the
        # second is the +20% revision. A stub that answered 1200 both times
        # would compare a run against itself and pass for the wrong reason.
        calls = {"n": 0}

        def t(url, payload, headers, timeout):
            captured["payload"] = payload
            calls["n"] += 1
            m = examples.build("plant_selection")
            if calls["n"] > 1:
                for c in m["constraints"]:
                    if c["name"] == "demand":
                        c["lower_bound"] = 1200.0
            return 200, envelope(json.dumps(m))

        a = make_app(client=make_client(t), key="k")
        first = a.route_formulate({"message": "three plants, demand 1000", "session": "s1"})
        self.assertTrue(first["ok"])
        run_a = solve_now(a, first["model"])

        second = a.route_formulate({"message": "what if demand increases by 20%?",
                                    "session": "s1"})
        self.assertTrue(second["ok"])
        sent = json.dumps(captured["payload"])
        self.assertIn("COMPLETE modified model", sent)
        self.assertIn("plant_selection", sent)

        run_b = solve_now(a, second["model"])
        self.assertEqual(run_b.result["status"], "OPTIMAL")
        self.assertTrue(run_b.result["checker_passed"])
        self.assertGreater(run_b.result["objective"], run_a.result["objective"])

    def test_solver_failure_after_a_successful_ai_call(self):
        """The AI succeeds, the model validates, and the model is genuinely
        infeasible. That must be reported as infeasible -- not as a failure
        of the AI, and not as a solution."""
        infeasible = json.loads(json.dumps(GOOD_MODEL))
        infeasible["constraints"][0]["lower_bound"] = 900  # demand above capacity
        t = transport_returning(json.dumps(infeasible))
        a = make_app(client=make_client(t), key="k")
        out = a.route_formulate({"message": "impossible demand"})
        self.assertTrue(out["ok"])
        run = solve_now(a, out["model"])
        self.assertEqual(run.state, "finished")
        self.assertIn(run.result["status"], ("INFEASIBLE", "ITERATION_LIMIT"))
        self.assertIsNone(run.result["objective"])


# ------------------------------------------------------------ offline demo

class OfflineDemo(unittest.TestCase):
    """Zero API calls, zero internet, zero LLM."""

    def test_every_example_validates_and_solves_verified(self):
        a = make_app(key="")  # its transport raises if anything tries to use it
        for item in examples.catalog():
            with self.subTest(example=item["id"]):
                out = a.route_example_model(item["id"])
                self.assertNotIsInstance(out, tuple, f"{item['id']} failed validation")
                run = solve_now(a, out["model"], time_limit=25.0)
                self.assertEqual(run.state, "finished", run.error)
                r = run.result
                self.assertEqual(r["status"], "OPTIMAL", item["id"])
                self.assertTrue(r["checker_passed"],
                                f"{item['id']}: independent checker did not pass")
                if r["is_mip"]:
                    self.assertTrue(r["proved_optimal"], f"{item['id']}: not proved optimal")

    def test_mip_emits_real_node_events(self):
        a = make_app(key="")
        out = a.route_example_model("plant_selection")
        run = solve_now(a, out["model"])
        evs = run.events()
        self.assertGreater(len(evs), 0)
        self.assertEqual(evs[0]["outcome"], "root")
        for e in evs:
            self.assertIn(e["outcome"], ("root", "branched", "integer", "infeasible",
                                          "dominated", "gap-cut", "unresolved"))
            self.assertGreaterEqual(e["depth"], 0)

    def test_lp_has_no_search_tree(self):
        a = make_app(key="")
        out = a.route_example_model("transportation")
        run = solve_now(a, out["model"])
        self.assertFalse(run.result["is_mip"])
        self.assertEqual(run.events(), [], "an LP must not report branch-and-bound nodes")

    def test_maximization_sign_is_restored(self):
        """A maximization is solved as a negated minimization. If the sign
        were not restored, the reported objective would be negative and the
        bound would be labelled the wrong way round."""
        a = make_app(key="")
        out = a.route_example_model("production_planning")
        run = solve_now(a, out["model"])
        self.assertEqual(run.result["objective_sense"], "maximize")
        self.assertGreater(run.result["objective"], 0)
        self.assertAlmostEqual(run.result["objective"], 9375.0, places=4)

    def test_edited_model_is_revalidated(self):
        a = make_app(key="")
        out = a.route_validate({"model": {"objective_sense": "minimize", "variables": [],
                                          "objective": [], "constraints": []}})
        self.assertFalse(out["valid"])
        self.assertEqual(len(a.runs), 0)


# ------------------------------------------------------------ round trip

class Serialization(unittest.TestCase):

    def test_no_non_finite_floats_survive_serialization(self):
        """Python's json emits a bare `Infinity` token for an IEEE infinity,
        which is not valid JSON and which the browser refuses to parse. The
        API sanitizes those to null; this asserts it, because the failure
        mode is a blank page rather than an error."""
        a = make_app(key="")
        out = a.route_example_model("production_planning")
        raw = json.dumps(app_module._sanitize(out), allow_nan=False)
        self.assertNotIn("Infinity", raw)
        self.assertNotIn("NaN", raw)

        def strict(c):
            raise AssertionError("bare constant " + c)
        json.loads(raw, parse_constant=strict)

    def test_model_survives_a_round_trip_through_the_api(self):
        """A model that came out of the API must be accepted back into it.
        This is what the UI does every time the user edits and re-runs."""
        a = make_app(key="")
        out = a.route_example_model("production_planning")
        wire = json.loads(json.dumps(app_module._sanitize(out), allow_nan=False))
        again = a.route_validate({"model": wire["model"]})
        self.assertTrue(again["valid"], again.get("errors"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
