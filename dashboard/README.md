# Inferno Console

A single-file interactive dashboard (`index.html`) — the "Solver Dashboard"
module from the team's SIH wireframe (`../../Inferno_SIH2026_Wireframe.pdf`,
module 6). Open it directly in a browser, or serve the folder
(`python3 -m http.server` from here).

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

No build step, no dependencies beyond two Google Fonts loaded over HTTPS.
Plain HTML/CSS/vanilla JS by design — this is a pitch/demo instrument, not
the start of a framework-based frontend; that decision can be revisited
once there's a real backend (Phase 2+) for it to talk to.
