#!/bin/bash
# ================================
# FILE      : build_petsc_sunbird2.sh
# PURPOSE   : ONE-OFF PETSc BUILD FOR CBS3D ON SUNBIRD2 (NO PETSc MODULE THERE).
#             OPTIMISED, SHARED, C/C++ ONLY, WITH HYPRE (BoomerAMG),
#             METIS/ParMETIS AND REFERENCE BLAS/LAPACK.
#
# USAGE     : RUN ON THE LOGIN NODE (configure DOWNLOADS PACKAGES):
#               mkdir -p ~/software && cd ~/software
#               nohup bash $CBS3D_CODE/jobs/sunbird/templates/build_petsc_sunbird2.sh \
#                   > build_petsc_sb2.log 2>&1 &
#               tail -f build_petsc_sb2.log
#             20-40 MINUTES. DONE WHEN THE LOG ENDS WITH "PETSC BUILD DONE".
#
# RESULT    : $HOME/software/petsc-sb2   (PETSC_DIR IN cbs3d_env.sh)
#             PETSc 3.26 + HYPRE 3.2 ON 29 OCT 2026.
# ================================

set -euo pipefail
type module > /dev/null 2>&1 || source /etc/profile
module purge
module load gnu/13.5.0-x86_64 mpi/openmpi/5.0.10-gnu-13.5.0-ucx-mt-x86_64 cmake/3.31.10
echo "mpicxx = $(command -v mpicxx)"
mpicxx --version | head -1

PREFIX="${PETSC_PREFIX:-$HOME/software/petsc-sb2}"
SRC="${PETSC_SRC:-$HOME/software/petsc-src}"
ARCH=arch-sb2-opt

[ -d "$SRC" ] || git clone --depth 1 -b release https://gitlab.com/petsc/petsc.git "$SRC"
cd "$SRC"

./configure --prefix="$PREFIX" PETSC_ARCH="$ARCH" \
    --with-cc=mpicc --with-cxx=mpicxx --with-fc=0 \
    --with-debugging=0 COPTFLAGS=-O3 CXXOPTFLAGS=-O3 --with-shared-libraries=1 \
    --download-f2cblaslapack --download-hypre --download-metis --download-parmetis
make PETSC_DIR="$PWD" PETSC_ARCH="$ARCH" MAKE_NP=8 all
make PETSC_DIR="$PWD" PETSC_ARCH="$ARCH" install
echo "PETSC BUILD DONE: $PREFIX"
