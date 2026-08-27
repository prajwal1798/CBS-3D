#!/bin/bash
# ==============================================================================
# Safe one-shot Sunbird gate for the wall-resolved NASA/TMR flat plate.
#
# 1. Build and run the deterministic Liu/Nithiarasu one-TET verifier.
# 2. Only if that passes, rebuild the MPI/PETSc production solver.
# 3. Submit a 100-iteration P40 smoke derived mechanically from the canonical
#    flatplate_sa_tmr_p40.slurm.  The 100k production script is not modified.
# ==============================================================================

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
VERIFY_BUILD=$CODE/build/liu-verification
PROD_EXE=$CODE/build/mpi-release/cbs3d_parallel
CANONICAL_JOB=$CODE/jobs/sunbird/flatplate_sa_tmr_p40.slurm
SMOKE_DIR=$CASE/jobs/generated

GCC_ROOT=/lustrehome/apps/compilers/gnu/12.1.0
MPI_ROOT=/apps/libraries/openmpi/4.1.6/el7/AVX512/gnu-12.1
export PETSC_DIR=/home/s.2337862/software/petsc-gcc12-install

cd "$CODE"

printf '===== SOURCE GATE =====\n'
BRANCH=$(git rev-parse --abbrev-ref HEAD)
HEAD=$(git rev-parse HEAD)
printf 'branch = %s\n' "$BRANCH"
printf 'HEAD   = %s\n' "$HEAD"
[[ "$BRANCH" == "integration/sa-mpi" ]] || {
    echo "FATAL: expected integration/sa-mpi" >&2
    exit 1
}
[[ -z "$(git status --porcelain)" ]] || {
    echo "FATAL: Sunbird worktree is dirty" >&2
    git status --short
    exit 1
}

printf '\n===== LIU/NITHIARASU DETERMINISTIC ELEMENT VERIFICATION =====\n'
module purge
module load cmake/3.31.10
export PATH="$GCC_ROOT/bin:$PATH"

rm -rf "$VERIFY_BUILD"
CXX="$GCC_ROOT/bin/g++" cmake \
    -S "$CODE/tests/liu" \
    -B "$VERIFY_BUILD" \
    -DCMAKE_BUILD_TYPE=Release
cmake --build "$VERIFY_BUILD" --parallel 4
ctest --test-dir "$VERIFY_BUILD" --output-on-failure

printf '\nLIU FORMULATION GATE: PASS\n'

printf '\n===== MPI/PETSC PRODUCTION BUILD =====\n'
module purge
module load cmake/3.31.10
export PATH="$MPI_ROOT/bin:$GCC_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$PETSC_DIR/lib:$MPI_ROOT/lib:$GCC_ROOT/lib64:${LD_LIBRARY_PATH:-}"
unset PETSC_ARCH
export OMPI_CC="$GCC_ROOT/bin/gcc"
export OMPI_CXX="$GCC_ROOT/bin/g++"

cmake --preset mpi-release
cmake --build --preset build-mpi-release --parallel 4
[[ -x "$PROD_EXE" ]] || {
    echo "FATAL: production executable missing after build: $PROD_EXE" >&2
    exit 1
}

printf '\n===== GENERATE 100-ITERATION P40 SMOKE =====\n'
[[ -f "$CANONICAL_JOB" ]] || {
    echo "FATAL: canonical TMR job missing: $CANONICAL_JOB" >&2
    exit 1
}
mkdir -p "$SMOKE_DIR" "$CASE/jobs/logs"
SMOKE_JOB="$SMOKE_DIR/flatplate_sa_liu_gate_${HEAD}.slurm"

python3 - "$CANONICAL_JOB" "$SMOKE_JOB" <<'PY'
from pathlib import Path
import sys

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
text = src.read_text(encoding="utf-8")

replacements = [
    ("#SBATCH --job-name=fp_sa_tmr_p40", "#SBATCH --job-name=fp_sa_gate"),
    ("#SBATCH --time=24:00:00", "#SBATCH --time=00:30:00"),
    ('time_values[0] = "100000"', 'time_values[0] = "100"'),
    ('time_values[4] = "1000"', 'time_values[4] = "100"'),
    ('lines[output_i] = "1  10  1000  0  1  10000  0.0  1  5000"',
     'lines[output_i] = "1  1  10  0  1  100  0.0  1  100"'),
    ("printf 'max iterations  = 100000\\n'", "printf 'max iterations  = 100\\n'"),
    ("printf 'SA stop gate    = 1.0e-7 after >=5000 iterations\\n'",
     "printf 'SA stop gate    = 1.0e-7 after >=100 iterations\\n'")
]

for old, new in replacements:
    count = text.count(old)
    if count != 1:
        raise SystemExit(
            "FATAL: canonical TMR job changed; expected exactly one occurrence of {!r}, found {}"
            .format(old, count))
    text = text.replace(old, new)

dst.write_text(text, encoding="utf-8")
print("Generated guarded smoke job:", dst)
PY

chmod 700 "$SMOKE_JOB"

printf '\n===== SUBMIT P40 SMOKE =====\n'
JOBID=$(sbatch --parsable "$SMOKE_JOB")
printf 'JOBID = %s\n' "$JOBID"
printf 'OUT   = %s\n' "$CASE/jobs/logs/sa_tmr_${JOBID}.out"
printf 'ERR   = %s\n' "$CASE/jobs/logs/sa_tmr_${JOBID}.err"
printf '\nLive output:\n'
printf 'tail --retry -F %s\n' "$CASE/jobs/logs/sa_tmr_${JOBID}.out"
printf '\nQueue:\n'
squeue -j "$JOBID" || true
