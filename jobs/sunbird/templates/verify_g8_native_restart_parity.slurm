#!/bin/bash -l
# G8a: compare native checkpoint versus legacy VTU restart.
# Replays one coupled, 2-rank iteration from the same G6 step-50000 state.
# Reuses EXACT production binary, MPI partition, .par and .matprop.
#SBATCH --job-name=cbs3d_g8a
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=3G
#SBATCH --time=00:10:00
#SBATCH --output=cbs3d_g8a_%j.out
#SBATCH --error=cbs3d_g8a_%j.err

set -euo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1
BASE="$CBS3D_SCRATCH/cases/uchannel_afc_audit"
ORIGINAL="$BASE/run_uchannel_29706"
G6ROOT="$CBS3D_SCRATCH/verification_g6/29763"
LEGACY_REPLAY="$G6ROOT/run_uchannel_29763/output/uchannel"
NATIVE="$LEGACY_REPLAY/restart/step_00050000"
MASTER="$G6ROOT/flow_one_step.par"
PAR="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.par"
MAT="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.matprop"
WORK="$CBS3D_SCRATCH/verification_g8a/${SLURM_JOB_ID:-manual}"
mkdir -p "$WORK"

for file in "$MASTER" "$PAR" "$MAT" \
    "$NATIVE/restart_manifest.txt" \
    "$NATIVE/uchannel_step_00050000_rank_0000.cbsrst" \
    "$NATIVE/uchannel_step_00050000_rank_0001.cbsrst"; do
    [[ -s "$file" ]] || { echo "G8a FATAL missing file: $file"; exit 1; }
done

hash_match() {
    local actual
    actual="$(sha256sum "$2" | awk '{print $1}')"
    echo "SHA256 $2 $actual"
    [[ "$1" == "$actual" ]] || { echo "G8a FATAL changed original file: $2"; exit 1; }
}
hash_match "96ca7134f6affc2cd17cdfb6b14ba70bbda9f243d7041ba196d70a00bc1c16fa" "$CBS3D_MPI_EXE"
hash_match "706b64a8e641a53294c393c23dba3c086cd104fd2961541af7cfe3d01839e6ac" "$PAR"
hash_match "928bc80463893a6d80c19eb9b5df27e5553ccd8b9821473935acf7e12641a3fb" "$MAT"
cmp -s "$MASTER" "$G6ROOT/run_uchannel_29763/partitions/rank_0000/uchannel_rank_0000.par" || {
    echo "G8a FATAL G6 reference .par differs from its staged copy"; exit 1; }

for name in \
    CBS3D_AVERAGE_START_ITERATION CBS3D_AVERAGE_START_TIME \
    CBS3D_CHECKPOINT_ROOT CBS3D_CONTINUATION_ITERATIONS \
    CBS3D_RESTART_VELOCITY_ARRAY CBS3D_SA_WALL_TREATMENT \
    CBS3D_CHT_WALL_TREATMENT \
    CBS3D_STEADY_ENERGY_MAXIT CBS3D_STEADY_ENERGY_RELAX \
    CBS3D_STEADY_ENERGY_TOL CBS3D_STEADY_ENERGY_RESIDUAL_TOL
do
    unset "$name"
done

export CBS3D_TEMPERATURE_AFC=0
export CBS3D_LTS_MONOTONE=1
export CBS3D_LTS_REFRESH_MARGIN=0.8
export CBS3D_STEADY_ENERGY=0
export CBS3D_CHECKPOINT_EVERY=0
export CBS3D_RESTART=1
export CBS3D_RESTART_FORMAT=native
export CBS3D_RESTART_ROOT="$NATIVE"
export CBS3D_RESTART_AS_INITIAL_CONDITION=0
export CBS3D_PAR="$MASTER"
export CBS3D_MATPROP="$MAT"

echo "DATE: $(date -Is)"
echo "SOURCE: $(git -C "$CBS3D_CODE" rev-parse HEAD)"
echo "Native restart: $NATIVE"
echo "Two-rank reference: $LEGACY_REPLAY"
echo "WORK: $WORK"
export SLURM_SUBMIT_DIR="$WORK"
bash "$CBS3D_CODE/jobs/sunbird/templates/run_mpi.slurm" \
    "$ORIGINAL/partitions" uchannel

CANDIDATE="$WORK/run_uchannel_${SLURM_JOB_ID:-manual}/output/uchannel"
python3 "$CBS3D_CODE/tools/compare_cbs_restart_parity_g8.py" \
    --legacy-vtu-dir "$LEGACY_REPLAY" \
    --native-vtu-dir "$CANDIDATE"
echo "G8A NATIVE RESTART REPLAY AND PARITY PASS"
