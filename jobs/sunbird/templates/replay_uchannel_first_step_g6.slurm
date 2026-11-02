#!/bin/bash -l
# G6 one-iteration replay from the original U-channel AFC-on step-50000.
# No original run/case/source inputs nor production solver code are modified.
#SBATCH --job-name=cbs3d_g6
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=3G
#SBATCH --time=00:10:00
#SBATCH --output=cbs3d_g6_%j.out
#SBATCH --error=cbs3d_g6_%j.err

set -euo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1

BASE="$CBS3D_SCRATCH/cases/uchannel_afc_audit"
ORIGINAL="$BASE/run_uchannel_29706"
SOURCE="$BASE/run_uchannel_29580/output/uchannel"
WORK="$CBS3D_SCRATCH/verification_g6/${SLURM_JOB_ID:-manual}"
mkdir -p "$WORK"
MASTER="$WORK/flow_one_step.par"
PAR="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.par"
MAT="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.matprop"

for file in "$PAR" "$MAT"; do
    [[ -s "$file" ]] || { echo "FATAL missing input: $file"; exit 1; }
done
[[ -d "$SOURCE" ]] || { echo "FATAL missing legacy restart directory $SOURCE"; exit 1; }

# Exact binary / rank0 case hashes recorded in Job 29706.
check_hash() {
    local actual
    actual="$(sha256sum "$2" | cut -d ' ' -f1)"
    echo "SHA256 $2 $actual"
    [[ "$actual" == "$1" ]] || { echo "FATAL hash mismatch: $2"; exit 1; }
}
check_hash "96ca7134f6affc2cd17cdfb6b14ba70bbda9f243d7041ba196d70a00bc1c16fa" "$CBS3D_MPI_EXE"
check_hash "706b64a8e641a53294c393c23dba3c086cd104fd2961541af7cfe3d01839e6ac" "$PAR"
check_hash "928bc80463893a6d80c19eb9b5df27e5553ccd8b9821473935acf7e12641a3fb" "$MAT"

# Duplicate the original .par; alter ONLY ntime and output cadence.
python3 - "$PAR" "$MASTER" <<'PY'
import pathlib
import sys
source = pathlib.Path(sys.argv[1])
destination = pathlib.Path(sys.argv[2])
lines = source.read_text().splitlines(keepends=True)
found = set()
for index, line in enumerate(lines[:-1]):
    head = line.strip().lower()
    if head.startswith("ntime transient_on"):
        values = lines[index+1].split()
        assert len(values) == 5, (head, values)
        assert values[0] == '50000', values
        assert values[1] == '0' and values[2] == '0', values
        values[0] = '1'
        lines[index+1] = '  '.join(values) + '\n'
        found.add('time')
    elif head.startswith("cbs3d output/monitor controls:"):
        values = lines[index+1].split()
        assert len(values) == 9, (head, values)
        values[0] = '1'
        values[1] = '1'
        values[2] = '1'
        values[4] = '1'
        values[5] = '1'
        lines[index+1] = '  '.join(values) + '\n'
        found.add('post')
assert found == {'time','post'}, found
destination.write_text(''.join(lines))
print("G6 .par changes: ntime=1, residual/console/VTU cadence=1")
PY

# Isolate this job from unrelated inherited options.
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
export CBS3D_RESTART_FORMAT=legacy_vtu
export CBS3D_RESTART_ROOT="$SOURCE"
export CBS3D_RESTART_ITERATION=50000
export CBS3D_RESTART_AS_INITIAL_CONDITION=0
export CBS3D_PAR="$MASTER"
export CBS3D_MATPROP="$MAT"

echo "G6 first-step replay: $(date -Is)"
echo "Git commit: $(git -C "$CBS3D_CODE" rev-parse HEAD)"
echo "Restart: $SOURCE step 50000; AFC=0; ilots=0; first iteration 50001"
echo "Original Job 29706: $ORIGINAL"
echo "Verification directory: $WORK"

# Existing trusted staging launcher puts a private copy of rank partitions
# under WORK; it never writes to ORIGINAL.
export SLURM_SUBMIT_DIR="$WORK"
bash "$CBS3D_CODE/jobs/sunbird/templates/run_mpi.slurm" \
    "$ORIGINAL/partitions" uchannel

REPLAY="$WORK/run_uchannel_${SLURM_JOB_ID:-manual}"
VTU="$REPLAY/output/uchannel"
echo "=== First-step output ==="
grep -E "Iteration [0-9]+/|Tmin=|stop reason|Run complete" "$REPLAY/solver.out" | tail -n 12 || true
echo "=== Cold-node comparison ==="
python3 "$CBS3D_CODE/tools/audit_temperature_undershoot_vtu.py" \
    --vtu-dir "$VTU" \
    --steps 50000,50001 \
    --reference-temperature 298.15 \
    --cold-tolerance 0.01

echo "G6 FIRST-STEP REPLAY COMPLETE: $VTU"
