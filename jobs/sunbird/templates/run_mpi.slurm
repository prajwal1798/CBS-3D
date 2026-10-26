#!/bin/bash -l
# ================================
# FILE      : run_mpi.slurm
# PURPOSE   : CBS3D DOMAIN-DECOMPOSED RUN WITH cbs3d_parallel.
#             ONE MPI RANK PER CORE, ONE RANK PER rank_XXXX FOLDER.
#
# USAGE     : sbatch run_mpi.slurm PARTITION_ROOT CASE_NAME
#             sbatch --nodes=2 --ntasks=80 run_mpi.slurm ...   (80 RANKS)
#
#   PARTITION_ROOT  FOLDER WITH rank_0000 ... rank_NNNN, MADE BY
#                   tools/export_rank_local_cbs.py --parts N
#                   (N MUST EQUAL --ntasks)
#
# OUTPUT    : <SUBMIT DIR>/run_<CASE>_<JOBID>/output/<CASE>/
#             (.pvtu + per-rank .vtu, <CASE>_distributed_residuals.csv)
#
# NOTES     : THE PARTITIONS ARE HARD-LINKED INTO THE RUN FOLDER, SO EDITING
#             A .par FOR ONE RUN NEVER CHANGES THE MASTER COPY.
# ================================
#SBATCH --job-name=cbs3d_mpi
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
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
for f in "$RUN_DIR"/partitions/rank_*/*.par; do
    cp --remove-destination "$(readlink -f "$f")" "$f.tmp" && mv "$f.tmp" "$f"
done

{
    echo "date      : $(date -Is)"
    echo "job       : $SLURM_JOB_ID, $SLURM_JOB_NUM_NODES node(s), $SLURM_NTASKS ranks"
    echo "partitions: $PART_ROOT"
    echo "exe       : $CBS3D_MPI_EXE"
    echo "git commit: $(git -C "$CBS3D_CODE" rev-parse HEAD 2>/dev/null || echo unknown)"
    sha256sum "$CBS3D_MPI_EXE" "$RUN_DIR/partitions/rank_0000/${CASE}_rank_0000.par"
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
