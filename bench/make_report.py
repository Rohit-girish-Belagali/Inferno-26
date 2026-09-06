#!/usr/bin/env python3
"""Turns bench/results.csv into the Phase 5 benchmark report.

BUILD_PLAN_V2.md's Phase 5 gate is "a stranger clones the repo and
reproduces our headline numbers using only the README", so this reads only
the committed results file and prints every number the README claims,
computed rather than transcribed. If a headline number in the docs and one
here ever disagree, this one is right.

bench/results.csv is committed for exactly that reason: it is the evidence
behind the claims, readable the moment the repo is cloned and without
waiting on a download and a full benchmark run. Regenerating it with
bench/run_netlib.py overwrites it, which is the point — a reader compares
their own run against ours rather than taking the committed numbers on
trust.

The performance profile is the Dolan-More construction the plan names. Its
point is that a single mean runtime hides the distribution: a solver can
win on average while failing badly somewhere. The profile plots, for each
factor t, the fraction of instances solved within t times the best known
time for that instance -- so the value at t=1 is "how often is this the
fastest" and the right-hand asymptote is "how often does it finish at all".
With one solver the reference time is its own best, which makes the curve a
plain runtime distribution; that is stated here rather than dressed up as a
comparison the data cannot support.
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results.csv")


def load(path):
    if not os.path.exists(path):
        sys.exit(f"no results at {path} — run bench/run_netlib.py first")
    with open(path) as f:
        return list(csv.DictReader(f))


def main():
    rows = load(RESULTS)
    total = len(rows)
    verified = [r for r in rows if r["status"] == "OPTIMAL_VERIFIED"]
    wrong = [r for r in rows if r["status"] == "OPTIMAL_CHECKER_FAILED"]
    by_status = {}
    for r in rows:
        by_status[r["status"]] = by_status.get(r["status"], 0) + 1

    print("=" * 68)
    print("INFERNO SOLVER — NETLIB BENCHMARK REPORT")
    print("=" * 68)
    print()
    print(f"Instances attempted:        {total}")
    print(f"Verified optimal:           {len(verified)}  ({100.0*len(verified)/total:.1f}%)")
    print(f"Checker-rejected claims:    {len(wrong)}"
          f"{'   <-- WRONG ANSWERS' if wrong else '   (none — no wrong answers)'}")
    print()
    print("Outcome breakdown:")
    for k in sorted(by_status, key=lambda k: -by_status[k]):
        print(f"  {k:<28} {by_status[k]}")
    print()

    times = sorted(float(r["time_s"]) for r in verified)
    if times:
        tot = sum(times)
        print("Solve time over verified instances:")
        print(f"  total                       {tot:.2f}s")
        print(f"  mean                        {tot/len(times):.3f}s")
        print(f"  median                      {times[len(times)//2]:.3f}s")
        print(f"  slowest                     {times[-1]:.3f}s")
        # The geometric mean is the standard summary for benchmark ratios:
        # unlike the arithmetic mean it is not dominated by the one slowest
        # instance, so it describes the typical case rather than the worst.
        gm = math.exp(sum(math.log(max(t, 1e-6)) for t in times) / len(times))
        print(f"  geometric mean              {gm:.4f}s")
        print()

    print("Slowest ten verified instances:")
    for r in sorted(verified, key=lambda r: -float(r["time_s"]))[:10]:
        print(f"  {r['name']:<14} {float(r['time_s']):>9.3f}s   obj = {r['objective']}")
    print()

    unsolved = [r for r in rows if r["status"] != "OPTIMAL_VERIFIED"]
    if unsolved:
        print("Not solved:")
        for r in sorted(unsolved, key=lambda r: r["name"]):
            print(f"  {r['name']:<14} {r['status']}")
        print()

    # --- Dolan-More performance profile ---
    print("Performance profile (Dolan-More)")
    print("  fraction of instances solved within a factor t of this solver's")
    print("  own best time. With a single solver the reference is its own")
    print("  best, so this is a runtime distribution, not a comparison.")
    print()
    if times:
        best = max(min(times), 1e-6)
        ratios = sorted(max(t, 1e-6) / best for t in times)
        print(f"  {'factor t':>12}   {'solved':>7}   profile")
        for t in [1, 2, 5, 10, 100, 1000, 10000, 100000]:
            k = sum(1 for r in ratios if r <= t)
            frac = k / total
            bar = "#" * int(round(frac * 40))
            print(f"  {t:>12}   {frac:>6.1%}   {bar}")
    print()
    print("Reproduce with:")
    print("  bash bench/download_netlib.sh")
    print("  python3 bench/run_netlib.py --solver revised --presolve")
    print("  python3 bench/make_report.py")
    print("  bash bench/verify_clean_room.sh")


if __name__ == "__main__":
    main()
