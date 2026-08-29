#!/bin/bash
# Build, verify and submit the fresh NASA/TMR SA-neg flat-plate P40 run.

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
EXE=$CODE/build/mpi-release/cbs3d_parallel
JOB=$CODE/jobs/sunbird/flatplate_sa_neg_tmr_p40.slurm

GCC_ROOT=/lustrehome/apps/compilers/gnu/12.1.0
MPI_ROOT=/apps/libraries/openmpi/4.1.6/el7/AVX512/gnu-12.1
PETSC_DIR=/home/s.2337862/software/petsc-gcc12-install

cd "$CODE"

[[ "$(git rev-parse --abbrev-ref HEAD)" == "integration/sa-mpi" ]] || {
    echo "FATAL: expected integration/sa-mpi" >&2
    exit 1
}

if ! grep -q 'NASA/TMR SA-neg, Allmaras-Johnson-Spalart' \
    src/assembly/SpalartAllmarasAssemblyLiu.cpp; then
    echo "FATAL: SA-neg production patch is not applied." >&2
    echo "Run: python3 tools/apply_sa_neg_production_branch.py" >&2
    exit 1
fi

if [[ -n "$(git status --porcelain)" ]]; then
    echo "FATAL: worktree is dirty. Commit the verified SA-neg patch before submitting." >&2
    git status --short
    exit 1
fi

module purge
module load cmake/3.31.10
export PATH="$MPI_ROOT/bin:$GCC_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$PETSC_DIR/lib:$MPI_ROOT/lib:$GCC_ROOT/lib64:${LD_LIBRARY_PATH:-}"
export PETSC_DIR
unset PETSC_ARCH
export OMPI_CC="$GCC_ROOT/bin/gcc"
export OMPI_CXX="$GCC_ROOT/bin/g++"

printf '===== SERIAL SA-neg VERIFICATION BUILD =====\n'
cmake -S . -B build/sa-neg-verify \
    -DCMAKE_BUILD_TYPE=Release \
    -DCBS3D_ENABLE_OPENMP=OFF \
    -DCBS3D_ENABLE_MPI=OFF \
    -DCBS3D_ENABLE_PETSC=OFF \
    -DCBS3D_BUILD_TESTS=ON \
    -DCBS3D_BUILD_SOURCE_VALIDATION=OFF
cmake --build build/sa-neg-verify --target sa_neg_element_verification --parallel 4
ctest --test-dir build/sa-neg-verify -R '^sa_neg_negative_branch$' --output-on-failure

printf '\n===== MPI/PETSc PRODUCTION BUILD =====\n'
cmake --preset mpi-release
cmake --build --preset build-mpi-release --parallel 4
[[ -x "$EXE" ]] || { echo "FATAL: executable missing: $EXE" >&2; exit 1; }

printf '\n===== SUBMIT FRESH SA-neg RUN =====\n'
mkdir -p "$CASE/jobs/logs"
JOBID=$(sbatch --parsable "$JOB")
echo "JOBID=$JOBID"
echo "RESULT_ROOT=$CASE/runs/sa_neg_tmr_wall_resolved/job_${JOBID}"
echo "LIVE:"
echo "tail --retry -F $CASE/runs/sa_neg_tmr_wall_resolved/job_${JOBID}/run/solver.out"
squeue -j "$JOBID"
