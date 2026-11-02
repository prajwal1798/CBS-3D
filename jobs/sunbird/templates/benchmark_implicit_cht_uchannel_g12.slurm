#!/bin/bash -l
# G12a: benchmark existing PETSc implicit steady CHT solver on the actual
# U-channel 2-rank step-50000 frozen velocity, Galerkin vs AFC.
# Standalone scratch build and run; never overwrites deployed executable.
#SBATCH --job-name=cbs3d_implicit_cht
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:40:00
#SBATCH --output=cbs3d_implicit_cht_%j.out
#SBATCH --error=cbs3d_implicit_cht_%j.err

set -euo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1
BASE="$CBS3D_SCRATCH/cases/uchannel_afc_audit"
ORIGINAL="$BASE/run_uchannel_29706"
REF="$BASE/run_uchannel_29580/output/uchannel"
PAR="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.par"
MAT="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.matprop"
WORK="$CBS3D_SCRATCH/benchmark_implicit_cht_g12/${SLURM_JOB_ID:-manual}"
BUILD="$WORK/build"
mkdir -p "$WORK"

for f in "$PAR" "$MAT" "$REF/uchannel_step_00050000_rank_0000.vtu" \
  "$REF/uchannel_step_00050000_rank_0001.vtu"; do
  [[ -s "$f" ]] || { echo "FATAL missing benchmark input: $f"; exit 1; }
done

# Independent candidate build, fixed source snapshot and provenance.
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
[[ -x "$EXE" ]] || { echo "FATAL benchmark build failed"; exit 1; }
echo "BENCHMARK GIT COMMIT: $(git -C "$CBS3D_CODE" rev-parse HEAD)"
sha256sum "$EXE" "$PAR" "$MAT"

# An isolated copy of the staging launcher uses the candidate executable.
python3 - "$CBS3D_CODE/jobs/sunbird/templates/run_mpi.slurm" "$WORK/launcher.sh" <<'PY'
from pathlib import Path
import sys
s=Path(sys.argv[1]).read_text()
x='source "$HOME/cbs3d_env.sh"'
assert s.count(x)==1
Path(sys.argv[2]).write_text(s.replace(x,x+'\nexport CBS3D_MPI_EXE="$CBS3D_G12_TEST_EXE"'))
PY

for name in \
    CBS3D_AVERAGE_START_ITERATION CBS3D_AVERAGE_START_TIME \
    CBS3D_CHECKPOINT_ROOT CBS3D_CONTINUATION_ITERATIONS \
    CBS3D_RESTART_VELOCITY_ARRAY CBS3D_SA_WALL_TREATMENT \
    CBS3D_CHT_WALL_TREATMENT CBS3D_AFC_SOLID \
    CBS3D_AFC_MIN_TEMPERATURE \
    CBS3D_STEADY_ENERGY_MAXIT CBS3D_STEADY_ENERGY_RELAX \
    CBS3D_STEADY_ENERGY_TOL CBS3D_STEADY_ENERGY_RESIDUAL_TOL \
    CBS3D_STEADY_ENERGY_ANDERSON_M CBS3D_STEADY_ENERGY_PSEUDOTIME_SHIFT \
    CBS3D_STEADY_ENERGY_ARMIJO_C1
do
    unset "$name"
done

export CBS3D_G12_TEST_EXE="$EXE"
export CBS3D_STEADY_ENERGY=1
export CBS3D_RESTART=1
export CBS3D_RESTART_FORMAT=legacy_vtu
export CBS3D_RESTART_ROOT="$REF"
export CBS3D_RESTART_ITERATION=50000
export CBS3D_RESTART_AS_INITIAL_CONDITION=0
export CBS3D_CHECKPOINT_EVERY=0
export CBS3D_PAR="$PAR"
export CBS3D_MATPROP="$MAT"
export CBS3D_LTS_MONOTONE=1
export CBS3D_LTS_REFRESH_MARGIN=0.8
export CBS3D_STEADY_ENERGY_TOL=1e-6
export CBS3D_STEADY_ENERGY_RESIDUAL_TOL=1e-6
export CBS3D_STEADY_ENERGY_MAXIT=150
export CBS3D_STEADY_ENERGY_RELAX=0.5
export CBS3D_STEADY_ENERGY_PSEUDOTIME_SHIFT=1.0
export CBS3D_AFC_LOG_EVERY=5000

failures=0
for variant in unlimited afc; do
    if [[ "$variant" == unlimited ]]; then
        export CBS3D_TEMPERATURE_AFC=0
    else
        export CBS3D_TEMPERATURE_AFC=1
    fi
    echo "=== G12a implicit $variant, AFC=$CBS3D_TEMPERATURE_AFC ==="
    mkdir -p "$WORK/$variant"
    export SLURM_SUBMIT_DIR="$WORK/$variant"
    export PETSC_OPTIONS="-energy_ksp_type gmres -energy_pc_type asm -energy_sub_pc_type ilu -energy_sub_pc_factor_levels 1 -energy_ksp_rtol 1e-10 -energy_ksp_max_it 4000"
    if bash "$WORK/launcher.sh" "$ORIGINAL/partitions" uchannel; then
        echo "G12a $variant: SOLVER FINISHED"
        RUN="$WORK/$variant/run_uchannel_${SLURM_JOB_ID:-manual}"
        grep -A 9 "Steady energy (implicit, frozen velocity)" "$RUN/solver.out" || true
        grep -E 'Run complete|stop reason|Tmin=|AFC.*iteration|AFC-PROBE' "$RUN/solver.out" | tail -n 25 || true
        echo "G12a $variant result: $RUN/output/uchannel"
    else
        echo "G12a $variant: SOLVER FAILED, inspect $WORK/$variant/run_uchannel_${SLURM_JOB_ID:-manual}/solver.err"
        failures=$((failures+1))
    fi
done
echo "G12a benchmark variants failed: $failures"
[[ "$failures" == 0 ]] || exit 1
echo "G12A IMPLICIT CHT BOTH-VARIANT SOLVE PASS"
