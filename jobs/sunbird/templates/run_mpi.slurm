#!/bin/bash -l
# ================================
# FILE      : run_mpi.slurm
# PURPOSE   : CBS3D DOMAIN-DECOMPOSED RUN WITH cbs3d_parallel.
#             ONE MPI RANK PER CORE, ONE RANK PER rank_XXXX FOLDER.
#
# USAGE     : sbatch run_mpi.slurm PARTITION_ROOT CASE_NAME
#             sbatch --nodes=2 --ntasks=256 run_mpi.slurm ...  (256 RANKS)
#             sbatch --ntasks=64 run_mpi.slurm ...             (64 RANKS, 1 NODE)
#             EASIER: bash submit_mpi.sh CASE_DIR CASE N   (PARTITIONS + RUN)
#
#   PARTITION_ROOT  FOLDER WITH rank_0000 ... rank_NNNN, MADE BY
#                   tools/export_rank_local_cbs.py --parts N
#                   (N MUST EQUAL --ntasks)
#
# OUTPUT    : <SUBMIT DIR>/run_<CASE>_<JOBID>/output/<CASE>/
#             (.pvtu + per-rank .vtu, <CASE>_distributed_residuals.csv)
#
# NOTES     : SUNBIRD2 NODE = 128 CORES (256 HYPERTHREADS). --hint=nomultithread
#             PUTS ONE RANK ON EACH PHYSICAL CORE.
#             THE PARTITIONS ARE HARD-LINKED INTO THE RUN FOLDER, SO EDITING
#             A .par FOR ONE RUN NEVER CHANGES THE MASTER COPY.
#             ONE .par / .matprop FOR ALL RANKS: export CBS3D_PAR=FILE AND/OR
#             CBS3D_MATPROP=FILE BEFORE sbatch; THE RUN COPIES IT INTO EVERY
#             RANK OF ITS PRIVATE COPY (THE PARTITION SET IS NOT CHANGED).
# ================================
#SBATCH --job-name=cbs3d_mpi
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --output=cbs3d_mpi_%j.out
#SBATCH --error=cbs3d_mpi_%j.err

set -euo pipefail

PART_ROOT="$(readlink -f "${1:?USAGE: sbatch run_mpi.slurm PARTITION_ROOT CASE_NAME}")"
CASE="${2:?USAGE: sbatch run_mpi.slurm PARTITION_ROOT CASE_NAME}"

source "$HOME/cbs3d_env.sh"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# ================================
# STEP 1: CHECK THE PARTITION MATCHES THE ALLOCATION
# ================================
NRANK=$(find "$PART_ROOT" -maxdepth 1 -type d -name 'rank_[0-9][0-9][0-9][0-9]' | wc -l)
[[ "$NRANK" -eq "$SLURM_NTASKS" ]] || {
    echo "FATAL: $NRANK rank folders but --ntasks=$SLURM_NTASKS" >&2; exit 1; }
for i in $(seq 0 $((NRANK - 1))); do
    r=$(printf '%04d' "$i")
    for ext in plt bco material par matprop mpi; do
        [[ -s "$PART_ROOT/rank_$r/${CASE}_rank_$r.$ext" ]] || {
            echo "FATAL: missing $PART_ROOT/rank_$r/${CASE}_rank_$r.$ext" >&2; exit 1; }
    done
done
[[ -x "$CBS3D_MPI_EXE" ]] || { echo "FATAL: missing $CBS3D_MPI_EXE" >&2; exit 1; }

# ================================
# STEP 2: STAGE A PRIVATE COPY (HARD LINKS: NO EXTRA DISK)
# ================================
RUN_DIR="$SLURM_SUBMIT_DIR/run_${CASE}_${SLURM_JOB_ID}"
mkdir -p "$RUN_DIR/partitions"
cp -al "$PART_ROOT/." "$RUN_DIR/partitions/"
# A .par EDITED LATER MUST NOT WRITE THROUGH THE HARD LINK: GIVE EACH RANK ITS OWN.
# CBS3D_PAR / CBS3D_MATPROP (OPTIONAL): ONE MASTER FILE COPIED TO EVERY RANK.
for ext in par matprop; do
    var="CBS3D_${ext^^}"
    master="${!var:-}"
    if [[ -n "$master" ]]; then
        [[ -s "$master" ]] || { echo "FATAL: $var=$master IS MISSING OR EMPTY" >&2; exit 1; }
        master="$(readlink -f "$master")"
        for d in "$RUN_DIR"/partitions/rank_*; do
            r="${d##*rank_}"
            cp --remove-destination "$master" "$d/${CASE}_rank_$r.$ext"
        done
        echo "MASTER $ext : $master -> ALL $NRANK RANKS"
    elif [[ "$ext" == par ]]; then
        for f in "$RUN_DIR"/partitions/rank_*/*.par; do
            cp --remove-destination "$(readlink -f "$f")" "$f.tmp" && mv "$f.tmp" "$f"
        done
    fi
done

{
    echo "date      : $(date -Is)"
    echo "job       : $SLURM_JOB_ID, $SLURM_JOB_NUM_NODES node(s), $SLURM_NTASKS ranks"
    echo "nodes     : $SLURM_JOB_NODELIST (tasks per node: ${SLURM_TASKS_PER_NODE:-?})"
    echo "partitions: $PART_ROOT"
    echo "exe       : $CBS3D_MPI_EXE"
    echo "checkpoint: every ${CBS3D_CHECKPOINT_EVERY:-0} iterations (0 = off; export CBS3D_CHECKPOINT_EVERY before sbatch)"
    echo "git commit: $(git -C "$CBS3D_CODE" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "master par: ${CBS3D_PAR:-none (partition copies)}"
    echo "master mat: ${CBS3D_MATPROP:-none (partition copies)}"
    echo "temp. AFC : ${CBS3D_TEMPERATURE_AFC:-0} (CBS3D_TEMPERATURE_AFC)"
    echo "LTS       : monotone ${CBS3D_LTS_MONOTONE:-0}, margin ${CBS3D_LTS_REFRESH_MARGIN:-0.8}"
    echo "average   : start time ${CBS3D_AVERAGE_START_TIME:-off}, start iteration ${CBS3D_AVERAGE_START_ITERATION:-off}"
    echo "restart   : ${CBS3D_RESTART:-0} ${CBS3D_RESTART_FORMAT:-native} ${CBS3D_RESTART_ROOT:-} ${CBS3D_RESTART_ITERATION:-} initial-condition=${CBS3D_RESTART_AS_INITIAL_CONDITION:-0}"
    sha256sum "$CBS3D_MPI_EXE" "$RUN_DIR/partitions/rank_0000/${CASE}_rank_0000.par" \
        "$RUN_DIR/partitions/rank_0000/${CASE}_rank_0000.matprop"
} > "$RUN_DIR/run_record.txt"

# ================================
# STEP 3: RUN
# ================================
cd "$RUN_DIR"
srun --mpi="$CBS3D_SRUN_MPI" --unbuffered \
    --ntasks="$SLURM_NTASKS" --cpus-per-task=1 \
    --distribution=block:block --cpu-bind=cores --kill-on-bad-exit=1 \
    stdbuf -oL -eL "$CBS3D_MPI_EXE" "$CASE" "$RUN_DIR/partitions" \
    > solver.out 2> solver.err

sacct -j "$SLURM_JOB_ID" --format=JobID,Elapsed,State,ExitCode,MaxRSS >> run_record.txt || true
grep -E "stop reason|Run complete" solver.out || true
echo "DONE: $RUN_DIR/output/$CASE"
