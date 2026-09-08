"""Runs the copilot server against a SCRIPTED stub of OpenRouter.

This is a test harness, not a product feature, and it exists for one
reason: without an API key the browser-side AI path -- formulation
rendering, clarification handling, the invalid-model banner, the
explanation bubble -- cannot be exercised at all. Substituting a scripted
transport lets every one of those be driven and looked at, while making it
completely obvious that no real language model was involved.

    python3 server/tests/fake_gateway_server.py --port 8421 [--mode MODE]

Modes: model (default), clarify, invalid, malformed, timeout, ratelimit.
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SERVER = os.path.dirname(HERE)
sys.path.insert(0, SERVER)

import app as app_module  # noqa: E402
import examples  # noqa: E402
import llm  # noqa: E402
from http.server import ThreadingHTTPServer  # noqa: E402


def build_transport(mode):
    # The what-if leg of the demo needs the stub to change its answer the
    # second time it is asked, or "increase demand 20%" returns the
    # original model and the comparison compares a run against itself.
    state = {"formulations": 0}

    def transport(url, payload, headers, timeout):
        body = json.dumps(payload)
        # Explanation requests carry the solver result; answer those with
        # prose that quotes the numbers it was actually given.
        if "inferno_result" in body:
            data = json.loads(body)
            blob = data["messages"][-1]["content"]
            start = blob.find("{")
            result = json.loads(blob[start:blob.rfind("}") + 1])["inferno_result"]
            obj = result.get("objective")
            proved = result.get("proved_optimal")
            picked = [v["name"] for v in result.get("variables", [])
                      if v["type"] == "binary" and abs(v["value"]) > 0.5]
            text = (
                f"Inferno {'proved this solution optimal' if proved else 'did not prove optimality'}.\n\n"
                f"**Objective: {obj}**\n\n"
                f"Selected: {', '.join(picked) if picked else 'no binary was chosen'}.\n\n"
                f"Nodes explored: {result.get('nodes')}. "
                f"Independent checker: {'PASS' if result.get('checker_passed') else 'FAIL'}.\n\n"
                "_(This text came from a scripted stub, not a real language model.)_")
            return 200, json.dumps({"choices": [{"message": {"content": text}}]}).encode()

        if mode == "timeout":
            raise TimeoutError("timed out")
        if mode == "ratelimit":
            return 429, b"{}"
        if mode == "malformed":
            content = "Sure! Here is the model: {oh no this is not json"
        elif mode == "clarify":
            content = json.dumps({"clarification":
                                  "What is the total demand you need to cover?"})
        elif mode == "invalid":
            bad = examples.build("plant_selection")
            bad["constraints"][0]["terms"][0]["variable"] = "variable_that_does_not_exist"
            bad["variables"][0]["upper_bound"] = 9  # a binary that is not binary
            content = json.dumps(bad)
        else:
            m = examples.build("plant_selection")
            state["formulations"] += 1
            last_user = ""
            for msg in payload.get("messages", []):
                if msg.get("role") == "user":
                    last_user = msg.get("content", "")
            if state["formulations"] > 1 and "20%" in last_user:
                for c in m["constraints"]:
                    if c["name"] == "demand":
                        c["lower_bound"] = round(c["lower_bound"] * 1.2, 6)
            content = json.dumps(m)
        return 200, json.dumps({"choices": [{"message": {"content": content}}]}).encode()

    return transport


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8421)
    ap.add_argument("--mode", default="model",
                    choices=["model", "clarify", "invalid", "malformed", "timeout", "ratelimit"])
    args = ap.parse_args()

    cfg = llm.Config({"OPENROUTER_API_KEY": "stub-key-not-real",
                      "OPENROUTER_TEXT_MODEL": "stub/scripted-gateway",
                      "OPENROUTER_MAX_ATTEMPTS": "2"})
    client = llm.OpenRouterClient(cfg, transport=build_transport(args.mode),
                                  sleep=lambda _s: None)
    application = app_module.App(config=cfg, client=client)
    client._log = application.log
    httpd = ThreadingHTTPServer(("127.0.0.1", args.port),
                                app_module.make_handler(application))
    print(f"SCRIPTED STUB GATEWAY — mode={args.mode} — http://127.0.0.1:{args.port}/")
    print("No real language model is contacted by this process.")
    httpd.serve_forever()


if __name__ == "__main__":
    main()
