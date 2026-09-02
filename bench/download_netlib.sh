#!/usr/bin/env bash
# Pulls the Netlib LP set from the authoritative source (netlib.org), which
# stores instances in a compact custom encoding, and decodes each one to
# plain-text MPS with the `emps` tool (see tools/NOTICE.md). Idempotent:
# re-running skips files already decoded.
#
# Usage: bench/download_netlib.sh [output_dir]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${1:-$SCRIPT_DIR/netlib}"
RAW_DIR="$OUT_DIR/raw"
MPS_DIR="$OUT_DIR/mps"
BASE_URL="https://www.netlib.org/lp/data"

mkdir -p "$RAW_DIR" "$MPS_DIR"

if [ ! -x "$SCRIPT_DIR/tools/emps" ]; then
  echo "building emps decompressor..."
  cc -O2 -o "$SCRIPT_DIR/tools/emps" "$SCRIPT_DIR/tools/emps.c"
fi
EMPS="$SCRIPT_DIR/tools/emps"

# The Netlib LP set's top-level file list, filtered to actual instances:
# excludes docs (changes, ascii, minos, nams.ps.gz), the emps tool itself,
# and the kennington/ subdirectory (a separate, non-canonical LP set).
#
# Also excludes mpc.src, stocfor3 and truss: on netlib.org these three are
# not emps-encoded MPS files but shar bundles containing a Fortran/C
# generator program that must be compiled and run to produce the instance.
# Wiring that up is future work (Phase 1.1 tracking: bench/tools/NOTICE.md);
# for now they are simply missing from the local set.
INSTANCES=(
  25fv47 80bau3b adlittle afiro agg agg2 agg3 bandm beaconfd blend
  bnl1 bnl2 boeing1 boeing2 bore3d brandy capri cycle czprob d2q06c
  d6cube degen2 degen3 dfl001 e226 etamacro fffff800 finnis fit1d fit1p
  fit2d fit2p forplan ganges gfrd-pnc greenbea greenbeb grow15 grow22
  grow7 israel kb2 lotfi maros maros-r7 modszk1 nesm perold
  pilot pilot.ja pilot.we pilot4 pilot87 pilotnov recipe sc105 sc205
  sc50a sc50b scagr25 scagr7 scfxm1 scfxm2 scfxm3 scorpion scrs8 scsd1
  scsd6 scsd8 sctap1 sctap2 sctap3 seba share1b share2b shell ship04l
  ship04s ship08l ship08s ship12l ship12s sierra stair standata standgub
  standmps stocfor1 stocfor2 tuff vtp.base wood1p woodw
)

ok=0
fail=0
for name in "${INSTANCES[@]}"; do
  mps_out="$MPS_DIR/${name}.mps"
  if [ -s "$mps_out" ]; then
    ok=$((ok + 1))
    continue
  fi
  raw_out="$RAW_DIR/${name}"
  if ! curl -sSL --fail --max-time 30 "$BASE_URL/$name" -o "$raw_out"; then
    echo "FETCH FAILED: $name"
    fail=$((fail + 1))
    continue
  fi
  if ! "$EMPS" "$raw_out" > "$mps_out" 2>/dev/null || [ ! -s "$mps_out" ]; then
    echo "DECODE FAILED: $name"
    rm -f "$mps_out"
    fail=$((fail + 1))
    continue
  fi
  ok=$((ok + 1))
done

echo "netlib set: $ok decoded, $fail failed, out of ${#INSTANCES[@]} — see $MPS_DIR"
