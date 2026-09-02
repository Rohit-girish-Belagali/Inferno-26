# Inferno Console + Bench Report

Two single-file pages, sharing one design system, from the team's SIH
wireframe (`../../Inferno_SIH2026_Wireframe.pdf`): `index.html` is module 6
("Solver Dashboard"), `report.html` is module 5 ("Benchmarking &
Performance Analytics"). Open either directly in a browser, or serve the
folder (`python3 -m http.server` from here) — they link to each other.

It walks the real Prepare → Solve → Verify pipeline from
`BUILD_PLAN_V2.md`'s system architecture diagram against a handful of
already-solved Netlib instances, using **real numbers** pulled from this
project's own bench runs (`bench/results.csv`) and the Phase 1.2 real-basis
verification (`la/manual_checks/`) — objective values, iteration counts,
solve times, and checker residuals are not invented.

The one deliberately fictional part: the GPU/PDLP lane in the concurrent
solve manager. Phase 3 (CUDA PDLP) hasn't been built — the dev machine has
no NVIDIA GPU, per `BUILD_PLAN_V2.md`'s own hardware-blocker note — so that
lane is a modeled projection, not a measurement, and is labeled "simulated"
everywhere it appears in the UI. Nothing here should be read as a
performance claim beyond what Phase 1.1/1.2 actually measured.

`report.html` embeds the complete real `bench/results.csv` (93/93
instances, no cherry-picking) as of the run that generated it — a sortable,
filterable, log-scale bar list plus a plain data table (kept alongside the
chart, not instead of it), with a per-row hover tooltip. Data goes stale
the moment the solver changes; regenerate it by rerunning
`bench/run_netlib.py`, converting the CSV to the embedded JSON array (see
git history for the conversion script), and republishing — there's no
live connection between the two.

No build step, no dependencies beyond two Google Fonts loaded over HTTPS.
Plain HTML/CSS/vanilla JS by design — these are pitch/demo instruments,
not the start of a framework-based frontend; that decision can be
revisited once there's a real backend (Phase 2+) for them to talk to.
