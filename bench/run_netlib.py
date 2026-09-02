#!/usr/bin/env python3
"""Runs the solver over every decoded Netlib instance and prints a score.

This is the harness named in BUILD_PLAN_V2.md Phase 1.1 and the standing
CI rule ("CI runs all Netlib instances on every commit and prints the
score"). It does not (yet) compare against published optimal objective
values — bench/netlib_optima.csv only has a couple of entries filled in
honestly rather than guessed, see that file's header. Until it is
populated, "PASS" here means: solver returned OPTIMAL and the independent
checker accepted the solution, nothing more.
"""
import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path

BENCH_DIR = Path(__file__).resolve().parent
REPO_ROOT = BENCH_DIR.parent
SOLVER = REPO_ROOT / "build" / "solver"
MPS_DIR = BENCH_DIR / "netlib" / "mps"
RESULTS_CSV = BENCH_DIR / "results.csv"
PER_INSTANCE_TIMEOUT_S = 60


def load_optima():
    path = BENCH_DIR / "netlib_optima.csv"
    optima = {}
    if not path.exists():
        return optima
    with open(path) as f:
        for row in csv.DictReader(f):
            if row["objective"]:
                optima[row["name"]] = float(row["objective"])
    return optima


def run_one(mps_path: Path, solver_name: str, use_presolve: bool):
    start = time.time()
    args = [str(SOLVER), str(mps_path), "--solver", solver_name]
    if use_presolve:
        args.append("--presolve")
    try:
        proc = subprocess.run(
            args, capture_output=True, text=True, timeout=PER_INSTANCE_TIMEOUT_S,
        )
    except subprocess.TimeoutExpired:
        return {"status": "TIMEOUT", "objective": None, "checker": "", "time": PER_INSTANCE_TIMEOUT_S}

    elapsed = time.time() - start
    status, objective, checker = "CRASH", None, ""
    for line in proc.stdout.splitlines():
        if line.startswith("checker:"):
            checker = line[len("checker:"):].strip()
        elif " status=" in line:
            for tok in line.split():
                if tok.startswith("status="):
                    status = tok[len("status="):]
                elif tok.startswith("objective="):
                    try:
                        objective = float(tok[len("objective="):])
                    except ValueError:
                        objective = None
    if proc.returncode != 0 and status == "CRASH":
        status = f"EXIT_{proc.returncode}"
    return {"status": status, "objective": objective, "checker": checker, "time": elapsed}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", choices=["revised", "dense"], default="revised",
                         help="which solve path to benchmark (default: revised, the real "
                              "Phase 2.1 sparse-LU simplex; dense is the Phase 1.1 throwaway)")
    parser.add_argument("--presolve", action="store_true",
                         help="run presolve (fixed-variable + empty-column removal) before "
                              "solving; only affects --solver revised")
    args = parser.parse_args()

    if not SOLVER.exists():
        print(f"solver binary not found at {SOLVER} — build it first "
              f"(cmake --build build, or see README.md)", file=sys.stderr)
        return 1
    instances = sorted(MPS_DIR.glob("*.mps"))
    if not instances:
        print(f"no .mps files in {MPS_DIR} — run bench/download_netlib.sh first", file=sys.stderr)
        return 1

    optima = load_optima()
    rows = []
    counts = {}
    verified = 0
    for mps_path in instances:
        name = mps_path.stem
        result = run_one(mps_path, args.solver, args.presolve)

        # Standing rule #2: the checker runs on every solve and its verdict
        # is what counts, not the solver's own claimed status. An OPTIMAL
        # whose checker did not PASS is a wrong answer, scored as such
        # rather than folded into "OPTIMAL" — that distinction is the whole
        # reason the checker exists.
        checker_passed = result["checker"].startswith("PASS")
        bucket = result["status"]
        if result["status"] == "OPTIMAL":
            bucket = "OPTIMAL_VERIFIED" if checker_passed else "OPTIMAL_CHECKER_FAILED"
        counts[bucket] = counts.get(bucket, 0) + 1
        if bucket == "OPTIMAL_VERIFIED":
            verified += 1

        match = ""
        if bucket == "OPTIMAL_VERIFIED" and name in optima and result["objective"] is not None:
            match = "yes" if abs(result["objective"] - optima[name]) <= 1e-6 * max(1.0, abs(optima[name])) else "NO"

        rows.append({
            "name": name, "status": bucket,
            "objective": result["objective"], "checker": result["checker"],
            "time_s": f"{result['time']:.3f}", "matches_published_optimum": match,
        })
        print(f"{name:16s} {bucket:22s} obj={result['objective']} "
              f"time={result['time']:.3f}s  {result['checker']}")

    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    total = len(instances)
    print(f"\nscore: {verified}/{total} OPTIMAL_VERIFIED  ({dict(counts)})")
    if counts.get("OPTIMAL_CHECKER_FAILED"):
        print(f"*** {counts['OPTIMAL_CHECKER_FAILED']} instance(s) solved to a claimed optimum "
              f"the checker REJECTED — see bench/results.csv, these are wrong answers, not noise ***")
    print(f"results written to {RESULTS_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
