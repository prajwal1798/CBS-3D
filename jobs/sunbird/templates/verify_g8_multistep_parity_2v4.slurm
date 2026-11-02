#!/bin/bash -l
# G8c: actual U-channel 100-step coupled solver 2-rank vs 4-rank parity.
# Same original binary and physical inputs; independent partitions;
# both start from the *same* step-50000 global fields.
# Do not alter production numerical source, original input meshes, or run files.
#SBATCH --job-name=cbs3d_g8c
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=3G
#SBATCH --time=00:20:00
#SBATCH --output=cbs3d_g8c_%j.out
#SBATCH --error=cbs3d_g8c_%j.err

set -euo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1
BASE="$CBS3D_SCRATCH/cases/uchannel_afc_audit"
G6ROOT="$CBS3D_SCRATCH/verification_g6/29763"
G8BROOT="$CBS3D_SCRATCH/verification_g8b/29788"
ORIGINAL="$BASE/run_uchannel_29706"
PART2="$ORIGINAL/partitions"
PART4="$G8BROOT/partitions_4"
RESTART2="$G6ROOT/run_uchannel_29763/output/uchannel"
RESTART4="$G8BROOT/restart_vtu_step_50000_four_rank"
PAR="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.par"
MAT="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.matprop"
WORK="$CBS3D_SCRATCH/verification_g8c/${SLURM_JOB_ID:-manual}"
MASTER="$WORK/flow_100_steps.par"
mkdir -p "$WORK/two" "$WORK/four"

hash_match() {
  local digest
  digest="$(sha256sum "$2" | awk '{print $1}')"
  echo "SHA256 $2 $digest"
  [[ "$1" == "$digest" ]] || { echo "G8c FATAL: original checksum changed: $2"; exit 1; }
}
hash_match "96ca7134f6affc2cd17cdfb6b14ba70bbda9f243d7041ba196d70a00bc1c16fa" "$CBS3D_MPI_EXE"
hash_match "706b64a8e641a53294c393c23dba3c086cd104fd2961541af7cfe3d01839e6ac" "$PAR"
hash_match "928bc80463893a6d80c19eb9b5df27e5553ccd8b9821473935acf7e12641a3fb" "$MAT"

for rank in 0 1; do
   r="$(printf '%04d' "$rank")"
   [[ -s "$RESTART2/uchannel_step_00050000_rank_$r.vtu" ]] || {
      echo "G8c FATAL missing two-rank source VTU"; exit 1; }
done
for rank in 0 1 2 3; do
   r="$(printf '%04d' "$rank")"
   [[ -s "$RESTART4/uchannel_step_00050000_rank_$r.vtu" ]] || {
      echo "G8c FATAL missing four-rank initial VTU"; exit 1; }
   [[ -s "$PART4/rank_$r/uchannel_rank_$r.mpi" ]] || {
      echo "G8c FATAL missing four-rank partition"; exit 1; }
done

# A fresh private parameter file; ONLY run length, monitor/output frequency
# and the minimum steady-iteration count are changed from original inputs.
python3 - "$PAR" "$MASTER" <<'PY'
from pathlib import Path
import sys
lines=Path(sys.argv[1]).read_text().splitlines()
have=set()
for i, line in enumerate(lines[:-1]):
    header=line.strip().lower()
    if header.startswith("ntime transient_on"):
        values=lines[i+1].split()
        assert len(values)==5 and values[1:3]==["0","0"], values
        values[0]="100"
        lines[i+1]="  ".join(values)
        have.add("ntime")
    elif header.startswith("cbs3d output/monitor controls:"):
        values=lines[i+1].split()
        assert len(values)==9, values
        values[0]="1"     # residual logging ON
        values[1]="100"   # log at 100 steps
        values[2]="100"   # console every 100
        values[4]="1"     # VTU ON
        values[5]="100"   # VTU initial and at step 50100
        values[6]="0"     # no fixed-time VTU
        values[8]="101"   # never stop early, one hundred iterations
        lines[i+1]="  ".join(values)
        have.add("output")
assert have=={"ntime","output"},have
Path(sys.argv[2]).write_text("\n".join(lines)+"\n")
print("G8c par: 100 iterations, VTU at 50000 and 50100, no early steady stop.")
PY

# Explicitly match AFC/restart/time schemes and unset unrelated overrides.
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
export CBS3D_RESTART_ITERATION=50000
export CBS3D_RESTART_AS_INITIAL_CONDITION=0
export CBS3D_PAR="$MASTER"
export CBS3D_MATPROP="$MAT"

echo "DATE: $(date -Is)"
echo "CODE COMMIT: $(git -C "$CBS3D_CODE" rev-parse HEAD)"
echo "G8c 2-v-4-rank, 100 coupled iterations; work: $WORK"

# Run sequentially on a single reserved 4-core job, using 2 then 4 slots.
# run_mpi.slurm honours SLURM_NTASKS and stages rank files privately.
(
   export SLURM_NTASKS=2
   export SLURM_SUBMIT_DIR="$WORK/two"
   export CBS3D_RESTART_ROOT="$RESTART2"
   bash "$CBS3D_CODE/jobs/sunbird/templates/run_mpi.slurm" "$PART2" uchannel
)
(
   export SLURM_NTASKS=4
   export SLURM_SUBMIT_DIR="$WORK/four"
   export CBS3D_RESTART_ROOT="$RESTART4"
   bash "$CBS3D_CODE/jobs/sunbird/templates/run_mpi.slurm" "$PART4" uchannel
)

VTU2="$WORK/two/run_uchannel_${SLURM_JOB_ID:-manual}/output/uchannel"
VTU4="$WORK/four/run_uchannel_${SLURM_JOB_ID:-manual}/output/uchannel"

python3 "$CBS3D_CODE/tools/verify_uchannel_partition_parity_g8b.py" compare \
   --reference "$VTU2" \
   --candidate "$VTU4" \
   --before 50000 --after 50100

echo "G8C 100-STEP 2-VS-4-RANK PARITY PASS"
