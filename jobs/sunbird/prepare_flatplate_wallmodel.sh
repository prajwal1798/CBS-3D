#!/bin/bash
# Synchronise the canonical Sunbird repository to integration/sa-mpi, rebuild
# the validated MPI/PETSc production executable, and audit the retained P40
# coarse flat-plate partition before submission.

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
PART=$CASE/partition/coarse/P0040/cbs
EXE=$CODE/build/mpi-release/cbs3d_parallel

GCC_ROOT=/lustrehome/apps/compilers/gnu/12.1.0
MPI_ROOT=/apps/libraries/openmpi/4.1.6/el7/AVX512/gnu-12.1
PETSC_DIR=/home/s.2337862/software/petsc-gcc12-install

printf '===== CANONICAL REPOSITORY =====\n'
cd "$CODE"
pwd
git status --short --branch
git branch --show-current
git rev-parse HEAD
git remote -v

if [[ -n "$(git status --porcelain)" ]]; then
    echo "ERROR: canonical Sunbird worktree is dirty; refusing to pull/build." >&2
    exit 1
fi

printf '\n===== SYNCHRONISE integration/sa-mpi =====\n'
git fetch origin integration/sa-mpi
git switch integration/sa-mpi
git pull --ff-only origin integration/sa-mpi

printf '\n===== UPDATED HEAD =====\n'
git rev-parse HEAD
git log -5 --oneline --decorate

printf '\n===== TOOLCHAIN =====\n'
module purge
module load cmake/3.31.10

export PATH="$MPI_ROOT/bin:$GCC_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$PETSC_DIR/lib:$MPI_ROOT/lib:$GCC_ROOT/lib64:${LD_LIBRARY_PATH:-}"
export PETSC_DIR
unset PETSC_ARCH

export OMPI_CC="$GCC_ROOT/bin/gcc"
export OMPI_CXX="$GCC_ROOT/bin/g++"

which cmake
which mpicxx
cmake --version | head -1
mpicxx --version | head -1
printf 'PETSC_DIR=%s\n' "$PETSC_DIR"

printf '\n===== CONFIGURE / BUILD MPI RELEASE =====\n'
cmake --preset mpi-release
cmake --build --preset build-mpi-release --parallel 4

[[ -x "$EXE" ]] || {
    echo "ERROR: expected executable not built: $EXE" >&2
    exit 1
}

printf '\n===== EXECUTABLE =====\n'
ls -lh "$EXE"

printf '\n===== FLAT-PLATE P40 PARTITION AUDIT =====\n'
[[ -d "$PART" ]] || {
    echo "ERROR: missing retained flat-plate partition root: $PART" >&2
    exit 1
}

rank_dirs=$(find "$PART" -maxdepth 1 -type d -name 'rank_[0-9][0-9][0-9][0-9]' | wc -l)
par_files=$(find "$PART" -mindepth 2 -maxdepth 2 -type f -name '*.par' | wc -l)
bco_files=$(find "$PART" -mindepth 2 -maxdepth 2 -type f -name '*.bco' | wc -l)
mpi_files=$(find "$PART" -mindepth 2 -maxdepth 2 -type f -name '*.mpi' | wc -l)

printf 'rank directories : %s\n' "$rank_dirs"
printf '.par files       : %s\n' "$par_files"
printf '.bco files       : %s\n' "$bco_files"
printf '.mpi files       : %s\n' "$mpi_files"

[[ "$rank_dirs" -eq 40 ]] || {
    echo "ERROR: expected 40 rank directories." >&2
    exit 1
}
[[ "$par_files" -eq 40 ]] || {
    echo "ERROR: expected 40 rank-local .par files." >&2
    exit 1
}
[[ "$bco_files" -eq 40 ]] || {
    echo "ERROR: expected 40 rank-local .bco files." >&2
    exit 1
}
[[ "$mpi_files" -eq 40 ]] || {
    echo "ERROR: expected 40 rank-local .mpi files." >&2
    exit 1
}

printf '\n===== REPRESENTATIVE PARAMETER / BC CONTROLS =====\n'
REP=$(find "$PART/rank_0000" -maxdepth 1 -type f -name '*.par' | head -1)
BCO=$(find "$PART/rank_0000" -maxdepth 1 -type f -name '*.bco' | head -1)
[[ -n "$REP" && -n "$BCO" ]] || {
    echo "ERROR: representative rank_0000 input files missing." >&2
    exit 1
}

grep -n -A1 -E 'ntime transient_on dtfixed dtfix iwrite|CBS3D timestep controls:|Spalart-Allmaras controls:|Spalart-Allmaras inlet/thermal|dimensional' "$REP" || true
printf '\nrank_0000 .bco:\n'
cat "$BCO"

printf '\n===== READY =====\n'
printf 'Executable : %s\n' "$EXE"
printf 'Partition  : %s\n' "$PART"
printf 'Submit     : sbatch %s/jobs/sunbird/flatplate_wallmodel_p40_startup.slurm\n' "$CODE"
