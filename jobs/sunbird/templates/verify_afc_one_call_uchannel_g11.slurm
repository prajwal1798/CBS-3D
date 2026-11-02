#!/bin/bash -l
# G11a: build corrected one-AFC-call binary in isolated scratch and
# perform one-step U-channel AFC-on experiment from G6 original restart.
# This does not overwrite the deployed executable or original case.
#SBATCH --job-name=cbs3d_afc_check
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:30:00
#SBATCH --output=cbs3d_afc_check_%j.out
#SBATCH --error=cbs3d_afc_check_%j.err

set -euo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1
BASE="$CBS3D_SCRATCH/cases/uchannel_afc_audit"
ORIGINAL="$BASE/run_uchannel_29706"
G6ROOT="$CBS3D_SCRATCH/verification_g6/29763"
REF="$G6ROOT/run_uchannel_29763/output/uchannel"
INPUT="$BASE/run_uchannel_29580/output/uchannel"
PAR="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.par"
MAT="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.matprop"
MASTER="$G6ROOT/flow_one_step.par"
WORK="$CBS3D_SCRATCH/verification_g11a/${SLURM_JOB_ID:-manual}"
BUILD="$WORK/build"
mkdir -p "$WORK"

for f in "$PAR" "$MAT" "$MASTER" \
  "$REF/uchannel_step_00050000_rank_0000.vtu" \
  "$REF/uchannel_step_00050001_rank_0000.vtu" \
  "$INPUT/uchannel_step_00050000_rank_0000.vtu"; do
  [[ -s "$f" ]] || { echo "FATAL missing $f"; exit 1; }
done
hash_eq() {
  local actual
  actual="$(sha256sum "$2" | awk '{print $1}')"
  [[ "$actual" == "$1" ]] || { echo "FATAL SHA mismatch $2"; exit 1; }
}
hash_eq 96ca7134f6affc2cd17cdfb6b14ba70bbda9f243d7041ba196d70a00bc1c16fa "$CBS3D_MPI_EXE"
hash_eq 706b64a8e641a53294c393c23dba3c086cd104fd2961541af7cfe3d01839e6ac "$PAR"
hash_eq 928bc80463893a6d80c19eb9b5df27e5553ccd8b9821473935acf7e12641a3fb "$MAT"

# A private build: original production binary remains untouched.
cmake -S "$CBS3D_CODE" -B "$BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCBS3D_ENABLE_MPI=ON \
  -DCBS3D_ENABLE_PETSC=ON \
  -DCBS3D_ENABLE_OPENMP=OFF \
  -DCBS3D_KEEP_LEGACY_EXECUTABLE_NAME=OFF \
  -DCBS3D_BUILD_SOURCE_VALIDATION=OFF \
  -DCBS3D_BUILD_TESTS=OFF

cmake --build "$BUILD" --target cbs3d_parallel -j 4
EXE="$BUILD/cbs3d_parallel"
[[ -x "$EXE" ]] || { echo "FATAL no compiled candidate binary"; exit 1; }
sha256sum "$EXE"
echo "SOURCE COMMIT $(git -C "$CBS3D_CODE" rev-parse HEAD)"

# run_mpi.slurm sources the user environment and resets CBS3D_MPI_EXE.
# Patch only a private copy of this launcher, NEVER the repository launcher.
python3 - "$CBS3D_CODE/jobs/sunbird/templates/run_mpi.slurm" "$WORK/launcher.sh" <<'PY'
import sys
from pathlib import Path
s=Path(sys.argv[1]).read_text()
match='source "$HOME/cbs3d_env.sh"'
assert s.count(match)==1, "launcher environment hook changed"
s=s.replace(match,match+'\nexport CBS3D_MPI_EXE="$CBS3D_G11_TEST_EXE"')
Path(sys.argv[2]).write_text(s)
PY

for name in \
  CBS3D_AVERAGE_START_ITERATION CBS3D_AVERAGE_START_TIME \
  CBS3D_CHECKPOINT_ROOT CBS3D_CONTINUATION_ITERATIONS \
  CBS3D_RESTART_VELOCITY_ARRAY CBS3D_SA_WALL_TREATMENT \
  CBS3D_CHT_WALL_TREATMENT CBS3D_AFC_SOLID \
  CBS3D_STEADY_ENERGY_MAXIT CBS3D_STEADY_ENERGY_RELAX \
  CBS3D_STEADY_ENERGY_TOL CBS3D_STEADY_ENERGY_RESIDUAL_TOL \
  CBS3D_AFC_MIN_TEMPERATURE
do unset "$name"; done

export CBS3D_G11_TEST_EXE="$EXE"
export CBS3D_TEMPERATURE_AFC=1
export CBS3D_AFC_LOG_EVERY=1
export CBS3D_STEADY_ENERGY=0
export CBS3D_CHECKPOINT_EVERY=0
export CBS3D_RESTART=1
export CBS3D_RESTART_FORMAT=legacy_vtu
export CBS3D_RESTART_ROOT="$INPUT"
export CBS3D_RESTART_ITERATION=50000
export CBS3D_RESTART_AS_INITIAL_CONDITION=0
export CBS3D_PAR="$MASTER"
export CBS3D_MATPROP="$MAT"
export CBS3D_LTS_MONOTONE=1
export CBS3D_LTS_REFRESH_MARGIN=0.8
export SLURM_SUBMIT_DIR="$WORK"
echo "G11a corrected AFC explicit first step from original hot restart"
bash "$WORK/launcher.sh" "$ORIGINAL/partitions" uchannel

SOLVER="$WORK/run_uchannel_${SLURM_JOB_ID:-manual}"
echo "=== AFC limiter invocation count ==="
grep "^Temperature AFC" "$SOLVER/solver.out" || true
count="$(grep -c "^Temperature AFC" "$SOLVER/solver.out" || true)"
[[ "$count" -eq 1 ]] || {
  echo "G11a FAIL: expected 1 limiter call, found $count"; exit 1; }

python3 "$CBS3D_CODE/tools/audit_afc_uchannel_first_step_g11.py" \
  --afc-off "$REF" \
  --afc-on "$SOLVER/output/uchannel"

echo "G11A CORRECTED AFC ONE-CALL AND BOUNDS PASS"
