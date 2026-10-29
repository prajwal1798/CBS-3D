#!/bin/bash
# ================================
# FILE      : build_cbs3d.sh
# PURPOSE   : BUILD THE TWO PRODUCTION SOLVERS WITH THE TOOLCHAIN IN
#             $HOME/cbs3d_env.sh:
#               build/serial-petsc-release/cbs3d_serial   (OpenMP + PETSc)
#               build/mpi-release/cbs3d_parallel          (MPI + PETSc)
#
# USAGE     : bash $CBS3D_CODE/jobs/sunbird/templates/build_cbs3d.sh [JOBS]
#             ABOUT 5 MINUTES WITH 8 JOBS; FINE ON THE LOGIN NODE.
#             REBUILD AFTER EVERY git pull.
# ================================

set -euo pipefail
source "$HOME/cbs3d_env.sh"
JOBS="${1:-8}"

cd "$CBS3D_CODE"
# THE PETSc PRESETS MUST USE THE SAME MPI WRAPPERS THAT BUILT PETSc.
export CC=mpicc CXX=mpicxx

for p in serial-petsc-release mpi-release; do
    if cmake --preset "$p" > "build_$p.log" 2>&1 &&
       cmake --build "build/$p" --parallel "$JOBS" >> "build_$p.log" 2>&1; then
        echo "OK   $p"
    else
        echo "FAIL $p (SEE $CBS3D_CODE/build_$p.log)"
        tail -30 "build_$p.log"
        exit 1
    fi
done

ls -lh "$CBS3D_SERIAL_EXE" "$CBS3D_MPI_EXE"
ldd "$CBS3D_MPI_EXE" | grep -E "petsc|libmpi|HYPRE" || true
echo "git commit: $(git rev-parse --short HEAD) ($(git rev-parse --abbrev-ref HEAD))"
