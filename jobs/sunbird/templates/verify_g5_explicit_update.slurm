#!/bin/bash -l
# G5b: first-step nonzero-velocity explicit thermal transport diagnostic.
# Runs the production EnergyAssembly.cpp operator in a small C++ CTest harness.
# Does not launch the full coupled CBS momentum/pressure solver.
#SBATCH --job-name=cbs3d_g5b
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:15:00
#SBATCH --output=cbs3d_g5b_%j.out
#SBATCH --error=cbs3d_g5b_%j.err

set -euo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1

WORK="$CBS3D_SCRATCH/verification_g5b/${SLURM_JOB_ID:-manual}"
BUILD="$WORK/build"
mkdir -p "$WORK"

echo "DATE: $(date -Is)"
echo "HOST: $(hostname)"
echo "SOURCE: $(git -C "$CBS3D_CODE" rev-parse HEAD)"
echo "CMAKE: $(cmake --version | head -n 1)"
echo "CXX: $(c++ --version | head -n 1)"
echo "BENCHMARK: nx=8/16, Pe=3/300, Courant_x=0.1/0.3/0.6; bounded hot perturbation"
echo "PRODUCTION NUMERICS MODIFIED: NO"

cmake -S "$CBS3D_CODE" -B "$BUILD" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCBS3D_ENABLE_MPI=OFF \
    -DCBS3D_ENABLE_PETSC=OFF \
    -DCBS3D_ENABLE_OPENMP=OFF \
    -DCBS3D_BUILD_SOURCE_VALIDATION=ON \
    -DCBS3D_BUILD_TESTS=OFF

cmake --build "$BUILD" \
    --target cbs3d_g5_explicit_update_audit -j 4

ctest --test-dir "$BUILD" \
    -R '^g5_explicit_update_monotonicity_audit$' \
    --output-on-failure -V

echo "G5B EXPLICIT UPDATE DIAGNOSTIC PASS"
