#!/bin/bash
# ================================
# FILE      : submit_mpi.sh
# PURPOSE   : ONE COMMAND FOR AN N-RANK CBS3D RUN:
#               1. IF CASE_DIR/partitions_N DOES NOT EXIST, SUBMIT
#                  export_partitions.slurm TO MAKE IT (A FEW MINUTES)
#               2. SUBMIT THE RUN, HELD UNTIL THE PARTITIONS PASS
#
# USAGE     : bash submit_mpi.sh CASE_DIR CASE N [MODE] [TIME]
#   CASE_DIR  FOLDER WITH serial/CASE.* (MADE BY cbs_case.sh / build_case.slurm)
#   N         MPI RANKS = PARTITIONS (ANY NUMBER; METIS SPLITS ANY N)
#   MODE      nodes  WHOLE EXCLUSIVE NODES, ceil(N/128) OF THEM (DEFAULT, FASTEST)
#             any    N CORES WHEREVER THEY ARE FREE (run_mpi_anycores.slurm)
#   TIME      SLURM TIME LIMIT (DEFAULT 04:00:00; compute ALLOWS 3-00:00:00)
#
# EXAMPLES  : bash $T/submit_mpi.sh $CBS3D_SCRATCH/cases/ldc3d ldc3d 256          # 2 NODES
#             bash $T/free_cores.sh                                               # SEE WHAT IS FREE
#             bash $T/submit_mpi.sh $CBS3D_SCRATCH/cases/ldc3d ldc3d 300 any      # SCATTERED CORES
# ================================

set -euo pipefail
[[ "$#" -ge 3 ]] || { sed -n '2,21p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

source "$HOME/cbs3d_env.sh"
CASE_DIR="$(readlink -f "$1")"
CASE="$2"
N="$3"
MODE="${4:-nodes}"
TIME="${5:-04:00:00}"
T="$CBS3D_CODE/jobs/sunbird/templates"
CORES_PER_NODE="${CBS3D_CORES_PER_NODE:-128}"
PART="partitions_$N"

[[ "$N" =~ ^[0-9]+$ && "$N" -ge 2 ]] || { echo "ERROR: N MUST BE AN INTEGER >= 2" >&2; exit 1; }
cd "$CASE_DIR"

# ================================
# STEP 1: PARTITIONS (REUSED WHEN COMPLETE)
# ================================
DEP=()
if [[ -d "$PART" ]]; then
    HAVE=$(find "$PART" -maxdepth 1 -type d -name 'rank_[0-9][0-9][0-9][0-9]' | wc -l)
    [[ "$HAVE" -eq "$N" ]] || {
        echo "ERROR: $CASE_DIR/$PART HAS $HAVE RANK FOLDERS, NOT $N (UNFINISHED?). DELETE IT AND RERUN." >&2
        exit 1; }
    echo "PARTITIONS : $CASE_DIR/$PART (EXISTING)"
else
    [[ -s "serial/$CASE.plt" ]] || { echo "ERROR: MISSING $CASE_DIR/serial/$CASE.plt" >&2; exit 1; }
    PJOB=$(sbatch --parsable --time=01:00:00 "$T/export_partitions.slurm" "serial/$CASE" "$N" "$PART")
    echo "PARTITIONS : JOB $PJOB -> $CASE_DIR/$PART"
    DEP=(--dependency="afterok:$PJOB")
fi

# ================================
# STEP 2: THE RUN
# ================================
case "$MODE" in
    nodes)
        NODES=$(( (N + CORES_PER_NODE - 1) / CORES_PER_NODE ))
        RJOB=$(sbatch --parsable "${DEP[@]}" --time="$TIME" \
            --nodes="$NODES" --ntasks="$N" "$T/run_mpi.slurm" "$PART" "$CASE")
        echo "RUN        : JOB $RJOB, $N RANKS ON $NODES WHOLE NODE(S)"
        ;;
    any)
        RJOB=$(sbatch --parsable "${DEP[@]}" --time="$TIME" \
            --ntasks="$N" "$T/run_mpi_anycores.slurm" "$PART" "$CASE")
        echo "RUN        : JOB $RJOB, $N RANKS ON ANY FREE CORES"
        ;;
    *)
        echo "ERROR: MODE MUST BE nodes OR any" >&2; exit 1 ;;
esac

echo "WATCH      : python3 $CBS3D_CODE/tools/watch_run.py $CASE_DIR/run_${CASE}_${RJOB}"
