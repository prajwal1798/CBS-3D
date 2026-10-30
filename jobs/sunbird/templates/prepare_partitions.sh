#!/bin/bash
# ================================
# FILE      : prepare_partitions.sh
# PURPOSE   : MAKE SEVERAL PARTITION SETS OF ONE SERIAL CASE IN ADVANCE
#             (partitions_40, partitions_80, ...), ONE SMALL SLURM JOB EACH,
#             SO A RUN CAN LATER START ON WHATEVER NUMBER OF CORES IS FREE.
#             SETS THAT ALREADY EXIST (WITH THE RIGHT NUMBER OF RANKS) ARE KEPT.
#
# USAGE     : bash prepare_partitions.sh CASE_DIR CASE N1 [N2 ...]
#   e.g.      bash $T/prepare_partitions.sh $CBS3D_SCRATCH/cases/ldc3d ldc3d 40 80 128 256 384
#
# NOTES     : EVERY SET GETS A COPY OF serial/CASE.par AT EXPORT TIME. EDIT
#             THE .par FIRST, OR COPY IT INTO EVERY rank_*/ AFTERWARDS.
# ================================

set -euo pipefail
[[ "$#" -ge 3 ]] || { sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

source "$HOME/cbs3d_env.sh"
CASE_DIR="$(readlink -f "$1")"
CASE="$2"
shift 2
T="$CBS3D_CODE/jobs/sunbird/templates"

cd "$CASE_DIR"
[[ -s "serial/$CASE.plt" ]] || { echo "ERROR: MISSING $CASE_DIR/serial/$CASE.plt" >&2; exit 1; }

for N in "$@"; do
    PART="partitions_$N"
    if [[ -d "$PART" ]]; then
        HAVE=$(find "$PART" -maxdepth 1 -type d -name 'rank_[0-9][0-9][0-9][0-9]' | wc -l)
        if [[ "$HAVE" -eq "$N" ]]; then
            echo "KEEP   : $PART ($HAVE ranks)"
        else
            echo "ERROR  : $PART HAS $HAVE RANK FOLDERS, NOT $N; DELETE IT AND RERUN" >&2
        fi
        continue
    fi
    JOB=$(sbatch --parsable --job-name="cbs3d_part$N" "$T/export_partitions.slurm" \
          "serial/$CASE" "$N" "$PART")
    echo "SUBMIT : $PART  (JOB $JOB, LOG cbs3d_part_$JOB.out)"
done
