#!/usr/bin/env bash
# Convert the completed P40 flat-plate VTU state to the native CBS3D restart
# format and submit a further 100k wall-resolved SA continuation.

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
SOURCE_JOB_ID=${1:-8197477}
SOURCE_JOB=$CASE/runs/sa_tmr_wall_resolved/job_${SOURCE_JOB_ID}
SOURCE_RUN=$SOURCE_JOB/run
SOURCE_OUTPUT=$SOURCE_RUN/output/flatplate
SOURCE_SOLVER_OUT=$SOURCE_RUN/solver.out
SOURCE_PART=$SOURCE_JOB/partitions
TOOL=$CODE/tools/flatplate_vtu_to_native_restart.py
SLURM=$CODE/jobs/sunbird/flatplate_sa_tmr_restart_p40.slurm
EXE=$CODE/build/mpi-release/cbs3d_parallel
ITERATION=100000
TAG=step_00100000
RESTART_DIR=$SOURCE_OUTPUT/restart_from_vtu/$TAG

[[ -x "$EXE" ]] || { echo "FATAL: solver executable missing: $EXE" >&2; exit 1; }
[[ -f "$TOOL" ]] || { echo "FATAL: restart converter missing: $TOOL" >&2; exit 1; }
[[ -f "$SLURM" ]] || { echo "FATAL: restart Slurm job missing: $SLURM" >&2; exit 1; }
[[ -d "$SOURCE_PART" ]] || { echo "FATAL: source partition missing: $SOURCE_PART" >&2; exit 1; }
[[ -d "$SOURCE_OUTPUT" ]] || { echo "FATAL: source output missing: $SOURCE_OUTPUT" >&2; exit 1; }
[[ -f "$SOURCE_SOLVER_OUT" ]] || { echo "FATAL: source solver.out missing: $SOURCE_SOLVER_OUT" >&2; exit 1; }

grep -q '^Run complete$' "$SOURCE_SOLVER_OUT" || {
    echo "FATAL: source job did not reach Run complete" >&2
    exit 1
}
grep -q 'Iteration 100000/100000' "$SOURCE_SOLVER_OUT" || {
    echo "FATAL: source job does not contain completed iteration 100000" >&2
    exit 1
}

piece_count=$(find "$SOURCE_OUTPUT" -maxdepth 1 -type f -name "flatplate_${TAG}_rank_*.vtu" | wc -l)
[[ "$piece_count" -eq 40 ]] || {
    echo "FATAL: expected 40 final VTU pieces at iteration 100000; found $piece_count" >&2
    exit 1
}

pvtu=$SOURCE_OUTPUT/flatplate_${TAG}.pvtu
[[ -f "$pvtu" ]] || { echo "FATAL: final PVTU missing: $pvtu" >&2; exit 1; }

# Recover the final reported LTS dt for restart metadata.  The continuation
# rebuilds the frozen LTS field from the restored state, but carrying the final
# value keeps the checkpoint metadata faithful to the completed source run.
DT=$(python3 - "$SOURCE_SOLVER_OUT" <<'PY'
from __future__ import print_function
import re
import sys

path = sys.argv[1]
pattern = re.compile(r"Iteration\s+100000/100000\s+dt=([^\s]+)")
value = None
with open(path, "r") as handle:
    for line in handle:
        match = pattern.search(line)
        if match:
            value = match.group(1)
if value is None:
    raise SystemExit("FATAL: could not recover final dt from source solver.out")
print(value)
PY
)

printf '===== SOURCE FLAT-PLATE STATE =====\n'
printf 'source job       : %s\n' "$SOURCE_JOB_ID"
printf 'source output    : %s\n' "$SOURCE_OUTPUT"
printf 'final iteration  : %s\n' "$ITERATION"
printf 'final dt         : %s\n' "$DT"
printf 'VTU pieces       : %s\n' "$piece_count"
printf 'wall model       : OFF in source benchmark\n'

printf '\n===== CREATE EXACT NATIVE RESTART =====\n'
python3 "$TOOL" \
    --input-dir "$SOURCE_OUTPUT" \
    --output-dir "$RESTART_DIR" \
    --case flatplate \
    --iteration "$ITERATION" \
    --mpi-size 40 \
    --physical-time 0.0 \
    --time-step "$DT" \
    --temperature-enabled 0 \
    --turbulence-enabled 1

[[ -f "$RESTART_DIR/restart_manifest.txt" ]] || {
    echo "FATAL: native restart manifest was not created" >&2
    exit 1
}
restart_count=$(find "$RESTART_DIR" -maxdepth 1 -type f -name 'flatplate_step_00100000_rank_*.cbsrst' | wc -l)
[[ "$restart_count" -eq 40 ]] || {
    echo "FATAL: native restart conversion produced $restart_count rank files" >&2
    exit 1
}

printf '\n===== NATIVE RESTART MANIFEST =====\n'
cat "$RESTART_DIR/restart_manifest.txt"

printf '\n===== SUBMIT ADDITIONAL 100K =====\n'
JOB_ID=$(sbatch \
    --parsable \
    --export=ALL,SOURCE_JOB_ID="$SOURCE_JOB_ID",RESTART_ROOT="$RESTART_DIR" \
    "$SLURM")

printf 'JOB ID       : %s\n' "$JOB_ID"
printf 'SOURCE JOB   : %s\n' "$SOURCE_JOB_ID"
printf 'RESTART ROOT : %s\n' "$RESTART_DIR"
printf 'OUTPUT       : %s/jobs/logs/sa_tmr_restart_%s.out\n' "$CASE" "$JOB_ID"
printf 'ERROR        : %s/jobs/logs/sa_tmr_restart_%s.err\n' "$CASE" "$JOB_ID"
printf 'RUN ROOT     : %s/runs/sa_tmr_wall_resolved/job_%s\n' "$CASE" "$JOB_ID"
printf '\nLive solver output after startup:\n'
printf 'tail -F %s/runs/sa_tmr_wall_resolved/job_%s/run/solver.out\n' "$CASE" "$JOB_ID"
printf '\nQueue:\n'
printf 'squeue -j %s\n' "$JOB_ID"
