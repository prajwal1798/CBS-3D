#!/bin/bash
# ==============================================================================
# Guarded Sunbird launcher for a 5000-iteration wall-resolved NASA/TMR SA run.
#
# This is deliberately NOT the 100k production benchmark.  It rebuilds the
# current MPI/PETSc solver, derives a temporary 5k job mechanically from the
# canonical P40 benchmark, preserves the canonical physics/BCs, writes VTU every
# 500 iterations, and submits one 40-rank single-node diagnostic job.
# ==============================================================================

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
CANONICAL_JOB=$CODE/jobs/sunbird/flatplate_sa_tmr_p40.slurm
GENERATED_DIR=$CASE/jobs/generated
PROD_EXE=$CODE/build/mpi-release/cbs3d_parallel

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

printf '\n===== LIU/NITHIARASU FORMULATION GATE =====\n'
module purge
module load cmake/3.31.10
export PATH="$GCC_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$GCC_ROOT/lib64:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="$GCC_ROOT/lib64:${LIBRARY_PATH:-}"
rm -rf "$CODE/build/liu-verification"
CXX="$GCC_ROOT/bin/g++" cmake \
    -S "$CODE/tests/liu" \
    -B "$CODE/build/liu-verification" \
    -DCMAKE_BUILD_TYPE=Release
cmake --build "$CODE/build/liu-verification" --parallel 4
ctest --test-dir "$CODE/build/liu-verification" --output-on-failure
printf 'LIU FORMULATION GATE: PASS\n'

printf '\n===== MPI/PETSC PRODUCTION BUILD =====\n'
module purge
module load cmake/3.31.10
export PATH="$MPI_ROOT/bin:$GCC_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$PETSC_DIR/lib:$MPI_ROOT/lib:$GCC_ROOT/lib64:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="$GCC_ROOT/lib64:${LIBRARY_PATH:-}"
unset PETSC_ARCH
export OMPI_CC="$GCC_ROOT/bin/gcc"
export OMPI_CXX="$GCC_ROOT/bin/g++"
cmake --preset mpi-release
cmake --build --preset build-mpi-release --parallel 4
[[ -x "$PROD_EXE" ]] || {
    echo "FATAL: production executable missing after build: $PROD_EXE" >&2
    exit 1
}

printf '\n===== GENERATE 5000-ITERATION P40 DIAGNOSTIC =====\n'
[[ -f "$CANONICAL_JOB" ]] || {
    echo "FATAL: canonical TMR job missing: $CANONICAL_JOB" >&2
    exit 1
}
mkdir -p "$GENERATED_DIR" "$CASE/jobs/logs"
JOB="$GENERATED_DIR/flatplate_sa_diag5000_${HEAD}.slurm"

python3 - "$CANONICAL_JOB" "$JOB" <<'PY'
from pathlib import Path
import sys

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
text = src.read_text(encoding='utf-8')

replacements = [
    ('#SBATCH --job-name=fp_sa_tmr_p40', '#SBATCH --job-name=fp_sa_5k'),
    ('#SBATCH --time=24:00:00', '#SBATCH --time=02:00:00'),
    ('time_values[0] = "100000"', 'time_values[0] = "5000"'),
    ('time_values[4] = "1000"', 'time_values[4] = "500"'),
    ('lines[output_i] = "1  10  1000  0  1  10000  0.0  1  5000"',
     'lines[output_i] = "1  10  100  0  1  500  0.0  1  5000"'),
    ("printf 'max iterations  = 100000\\n'", "printf 'max iterations  = 5000\\n'")
]

for old, new in replacements:
    count = text.count(old)
    if count != 1:
        raise SystemExit(
            'FATAL: canonical job changed; expected exactly one occurrence of {!r}, found {}'.format(old, count)
        )
    text = text.replace(old, new)

# The diagnostic must remain single-node and wall-resolved.
required = [
    '#SBATCH --nodes=1',
    '#SBATCH --ntasks=40',
    'export OMPI_MCA_btl="^openib"',
    'unset CBS3D_SA_WALL_TREATMENT',
    'time_values[0] = "5000"',
    'lines[output_i] = "1  10  100  0  1  500  0.0  1  5000"'
]
for token in required:
    if token not in text:
        raise SystemExit('FATAL: generated diagnostic lost required guard: ' + token)

# No accidental production-length run is allowed from this launcher.
if 'time_values[0] = "100000"' in text:
    raise SystemExit('FATAL: 100k iteration setting leaked into diagnostic job')

dst.write_text(text, encoding='utf-8')
print('Generated guarded 5k job:', dst)
PY

chmod 700 "$JOB"

printf '\n===== SUBMIT 5000-ITERATION P40 DIAGNOSTIC =====\n'
JOBID=$(sbatch --parsable "$JOB")
printf 'JOBID = %s\n' "$JOBID"
printf 'OUT   = %s\n' "$CASE/jobs/logs/sa_tmr_${JOBID}.out"
printf 'ERR   = %s\n' "$CASE/jobs/logs/sa_tmr_${JOBID}.err"
printf 'RUN   = %s\n' "$CASE/runs/sa_tmr_wall_resolved/job_${JOBID}"
printf '\nLive output:\n'
printf 'tail --retry -F %s\n' "$CASE/jobs/logs/sa_tmr_${JOBID}.out"
printf '\nQueue:\n'
squeue -j "$JOBID" || true
