#!/bin/bash
#SBATCH --job-name=cbs3d_build_matrix
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=cbs3d_build_matrix_%j.out
#SBATCH --error=cbs3d_build_matrix_%j.err

set -euo pipefail

export GCC_ROOT=/lustrehome/apps/compilers/gnu/12.1.0
export MPI_ROOT=/apps/libraries/openmpi/4.1.6/el7/AVX512/gnu-12.1
export PETSC_DIR=/home/s.2337862/software/petsc-gcc12-install
unset PETSC_ARCH

export PATH="${MPI_ROOT}/bin:${GCC_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${PETSC_DIR}/lib:${MPI_ROOT}/lib:${GCC_ROOT}/lib64:${LD_LIBRARY_PATH:-}"

module load cmake/3.31.10

SOURCE="${CBS3D_SOURCE:-${SLURM_SUBMIT_DIR}}"

if [[ ! -f "${SOURCE}/CMakeLists.txt" ]]; then
    echo "ERROR: CBS3D source directory is invalid: ${SOURCE}" >&2
    echo "Submit from the repository root or export CBS3D_SOURCE." >&2
    exit 1
fi

cd "${SOURCE}"

export CC="${GCC_ROOT}/bin/gcc"
export CXX="${GCC_ROOT}/bin/g++"

build_preset()
{
    local configure_preset="$1"
    local build_preset="$2"
    local executable="$3"
    local build_directory="${SOURCE}/build/${configure_preset}"

    echo
    echo "============================================================"
    echo "CONFIGURE: ${configure_preset}"
    echo "============================================================"

    cmake --preset "${configure_preset}"
    cmake --build --preset "${build_preset}" \
        --parallel "${SLURM_CPUS_PER_TASK}"

    # Production presets may build only the selected solver application.
    # Build and execute the shared Step-4 volumetric-source validation
    # explicitly so every Sunbird variant is checked against the same source
    # kernel before it is accepted as ready for use.
    cmake --build "${build_directory}" \
        --target cbs3d_source_validation \
        --parallel "${SLURM_CPUS_PER_TASK}"

    ctest --test-dir "${build_directory}" \
        --output-on-failure

    if [[ ! -x "${executable}" ]]; then
        echo "ERROR: expected executable was not produced: ${executable}" >&2
        exit 1
    fi

    ls -lh "${executable}"
}

build_preset \
    serial-release \
    build-serial-release \
    "${SOURCE}/build/serial-release/cbs3d_serial"

build_preset \
    serial-openmp-release \
    build-serial-openmp-release \
    "${SOURCE}/build/serial-openmp-release/cbs3d_serial"

# PETSc on Sunbird was built against the validated OpenMPI/GCC toolchain.
export CC="${MPI_ROOT}/bin/mpicc"
export CXX="${MPI_ROOT}/bin/mpicxx"

build_preset \
    serial-petsc-release \
    build-serial-petsc-release \
    "${SOURCE}/build/serial-petsc-release/cbs3d_serial"

build_preset \
    mpi-release \
    build-mpi-release \
    "${SOURCE}/build/mpi-release/cbs3d_parallel"

build_preset \
    hybrid-release \
    build-hybrid-release \
    "${SOURCE}/build/hybrid-release/cbs3d_parallel"

echo
echo "============================================================"
echo "LINKAGE CHECK: MPI RELEASE"
echo "============================================================"

ldd "${SOURCE}/build/mpi-release/cbs3d_parallel" | \
    grep -E 'petsc|mpi|open-rte|open-pal|gomp' || true

echo
echo "CBS3D PRODUCTION BUILD MATRIX COMPLETED"
