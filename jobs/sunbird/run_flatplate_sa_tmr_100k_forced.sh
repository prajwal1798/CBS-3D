#!/bin/bash
set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
CANONICAL=$CODE/jobs/sunbird/flatplate_sa_tmr_p40.slurm
GENERATED_DIR=$CASE/jobs/generated
EXE=$CODE/build/mpi-release/cbs3d_parallel

GCC_ROOT=/lustrehome/apps/compilers/gnu/12.1.0
MPI_ROOT=/apps/libraries/openmpi/4.1.6/el7/AVX512/gnu-12.1
export PETSC_DIR=/home/s.2337862/software/petsc-gcc12-install

cd "$CODE"

echo "===== SOURCE GATE ====="
BRANCH=$(git rev-parse --abbrev-ref HEAD)
HEAD=$(git rev-parse HEAD)
echo "branch = $BRANCH"
echo "HEAD   = $HEAD"
[[ "$BRANCH" == "integration/sa-mpi" ]] || { echo "FATAL: expected integration/sa-mpi" >&2; exit 1; }
[[ -z "$(git status --porcelain)" ]] || { echo "FATAL: Sunbird worktree is dirty" >&2; git status --short; exit 1; }

[[ -f "$CANONICAL" ]] || { echo "FATAL: canonical job missing: $CANONICAL" >&2; exit 1; }
mkdir -p "$GENERATED_DIR" "$CASE/jobs/logs"
JOB="$GENERATED_DIR/flatplate_sa_tmr_100k_forced_${HEAD}.slurm"

python3 - "$CANONICAL" "$JOB" <<'PY'
from pathlib import Path
import sys

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
text = src.read_text(encoding='utf-8')

old_name = '#SBATCH --job-name=fp_sa_tmr_p40'
new_name = '#SBATCH --job-name=fp_sa_100k'
if text.count(old_name) != 1:
    raise SystemExit('FATAL: canonical job-name anchor changed')
text = text.replace(old_name, new_name)

old_output = 'lines[output_i] = "1  10  1000  0  1  10000  0.0  1  5000"'
new_output = 'lines[output_i] = "1  10  1000  0  1  10000  0.0  1  100000"'
if text.count(old_output) != 1:
    raise SystemExit('FATAL: canonical steady-minimum anchor changed')
text = text.replace(old_output, new_output)

old_print = "printf 'SA stop gate    = 1.0e-7 after >=5000 iterations\\n'"
new_print = "printf 'SA stop gate    = disabled before iteration 100000\\n'"
if text.count(old_print) != 1:
    raise SystemExit('FATAL: canonical stop-gate print anchor changed')
text = text.replace(old_print, new_print)

required = [
    '#SBATCH --nodes=1',
    '#SBATCH --ntasks=40',
    '#SBATCH --ntasks-per-node=40',
    'time_values[0] = "100000"',
    'lines[output_i] = "1  10  1000  0  1  10000  0.0  1  100000"',
    'export OMPI_MCA_btl="^openib"',
    'unset CBS3D_SA_WALL_TREATMENT',
    'unset CBS3D_CHT_WALL_TREATMENT',
    'sa_values[:3] != ["1", "0", "0"]',
    'max iterations  = 100000',
]
for token in required:
    if token not in text:
        raise SystemExit('FATAL: forced-100k job lost required guard: ' + token)

if 'lines[output_i] = "1  10  1000  0  1  10000  0.0  1  5000"' in text:
    raise SystemExit('FATAL: early steady-stop minimum remains in forced-100k job')

dst.write_text(text, encoding='utf-8')
print('Generated exact 100000-iteration job:', dst)
PY

chmod 700 "$JOB"

echo "===== BUILD MPI/PETSC RELEASE ====="
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
[[ -x "$EXE" ]] || { echo "FATAL: executable missing after build: $EXE" >&2; exit 1; }

echo "===== FORCED 100K JOB AUDIT ====="
grep -F '#SBATCH --job-name=fp_sa_100k' "$JOB"
grep -F 'time_values[0] = "100000"' "$JOB"
grep -F 'lines[output_i] = "1  10  1000  0  1  10000  0.0  1  100000"' "$JOB"
grep -F 'unset CBS3D_SA_WALL_TREATMENT' "$JOB"
grep -F 'export OMPI_MCA_btl="^openib"' "$JOB"
echo "FORCED 100K JOB AUDIT: PASS"

echo "===== SUBMIT EXACT 100000-ITERATION P40 RUN ====="
JOBID=$(sbatch --parsable "$JOB")
echo "JOBID = $JOBID"
echo "OUT   = $CASE/jobs/logs/sa_tmr_${JOBID}.out"
echo "ERR   = $CASE/jobs/logs/sa_tmr_${JOBID}.err"
echo "RUN   = $CASE/runs/sa_tmr_wall_resolved/job_${JOBID}"
echo
echo "Live output:"
echo "tail --retry -F $CASE/jobs/logs/sa_tmr_${JOBID}.out"
echo
echo "Queue:"
squeue -j "$JOBID" || true
