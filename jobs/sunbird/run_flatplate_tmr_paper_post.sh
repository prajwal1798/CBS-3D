#!/usr/bin/env bash
set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
JOB=${1:-}

if [[ -z "$JOB" ]]; then
    JOB=$(find "$CASE/runs/sa_tmr_wall_resolved" -maxdepth 1 -type d -name 'job_*' -printf '%f\n' 2>/dev/null \
        | sed 's/^job_//' \
        | sort -n \
        | while read -r j; do
            [[ -f "$CASE/runs/sa_tmr_wall_resolved/job_${j}/run/output/flatplate/flatplate_step_00200000.pvtu" ]] && echo "$j"
          done \
        | tail -1)
fi

[[ -n "$JOB" ]] || { echo "FATAL: no completed 200000-iteration flat-plate continuation found" >&2; exit 1; }

JOB_ROOT=$CASE/runs/sa_tmr_wall_resolved/job_${JOB}
PART=$JOB_ROOT/partitions
OUT=$JOB_ROOT/run/output/flatplate
POST=$OUT/post_flatplate_sa_global
PLOTS=$POST/paper_plots
PRE=$CODE/tools/postprocess_flatplate_sa_pvtu_py36.py
PAPER=$CODE/tools/plot_flatplate_paper_tmr_py36.py

[[ -d "$PART" ]] || { echo "FATAL: partition root missing: $PART" >&2; exit 1; }
[[ -d "$OUT" ]] || { echo "FATAL: output directory missing: $OUT" >&2; exit 1; }
[[ -f "$PRE" ]] || { echo "FATAL: distributed postprocessor missing: $PRE" >&2; exit 1; }
[[ -f "$PAPER" ]] || { echo "FATAL: paper plotter missing: $PAPER" >&2; exit 1; }

printf '===== FLAT-PLATE TMR PAPER POSTPROCESS =====\n'
printf 'job      : %s\n' "$JOB"
printf 'partition: %s\n' "$PART"
printf 'output   : %s\n' "$OUT"
printf 'post-dir : %s\n' "$POST"

rm -rf "$POST"
mkdir -p "$POST"

python3 "$PRE" \
    --partition-root "$PART" \
    --output-dir "$OUT" \
    --post-dir "$POST" \
    --mpi-size 40 \
    --step 200000 \
    --rho 1.0 \
    --mu 2.0e-7 \
    --u-inf 1.0 \
    --wall-bc 530 \
    --cf-bins 120 \
    --stations 0.97008 1.90334

python3 "$PAPER" \
    --post-dir "$POST" \
    --uinf 1.0

printf '\n===== VALIDATION SUMMARY =====\n'
cat "$PLOTS/validation_summary.txt"

printf '\n===== GENERATED FIGURES =====\n'
find "$PLOTS" -maxdepth 1 -type f \( -name '*.png' -o -name '*.pdf' \) -printf '%f\n' | sort

printf '\nFLAT-PLATE TMR PAPER POSTPROCESS: PASS\n'
printf 'plots   : %s\n' "$PLOTS"
printf 'summary : %s/validation_summary.txt\n' "$PLOTS"
