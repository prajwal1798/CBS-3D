#!/bin/bash -l
#SBATCH --job-name=fp_sa_r100k
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/s.2337862/CBS3D/cases/flatplate/jobs/logs/sa_tmr_restart_%j.out
#SBATCH --error=/scratch/s.2337862/CBS3D/cases/flatplate/jobs/logs/sa_tmr_restart_%j.err

# Continue the wall-resolved NASA/TMR flat-plate SA benchmark for up to another
# 100000 iterations from an exact native checkpoint at global iteration 100000.
# The source job id and checkpoint directory are supplied by the launcher.

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
EXE=$CODE/build/mpi-release/cbs3d_parallel
ANALYZER=$CODE/examples/Flat_Plate_Turbulent/analyze_distributed_tmr.py

: "${SOURCE_JOB_ID:?SOURCE_JOB_ID must be exported by the launcher}"
: "${RESTART_ROOT:?RESTART_ROOT must be exported by the launcher}"

SOURCE_JOB=$CASE/runs/sa_tmr_wall_resolved/job_${SOURCE_JOB_ID}
SOURCE_PART=$SOURCE_JOB/partitions
RUN_ROOT=$CASE/runs/sa_tmr_wall_resolved
JOB_ROOT=$RUN_ROOT/job_${SLURM_JOB_ID}
STAGED=$JOB_ROOT/partitions
RUN_DIR=$JOB_ROOT/run
CHECKPOINT_ROOT=$RUN_DIR/output/flatplate/restart

GCC_ROOT=/lustrehome/apps/compilers/gnu/12.1.0
MPI_ROOT=/apps/libraries/openmpi/4.1.6/el7/AVX512/gnu-12.1
PETSC_DIR=/home/s.2337862/software/petsc-gcc12-install

module purge
export PATH="$MPI_ROOT/bin:$GCC_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$PETSC_DIR/lib:$MPI_ROOT/lib:$GCC_ROOT/lib64:${LD_LIBRARY_PATH:-}"
export PETSC_DIR
unset PETSC_ARCH
export OMPI_CC="$GCC_ROOT/bin/gcc"
export OMPI_CXX="$GCC_ROOT/bin/g++"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Wall-resolved benchmark only.  A wall law here would invalidate the test.
unset CBS3D_SA_WALL_TREATMENT
unset CBS3D_CHT_WALL_TREATMENT
unset CBS3D_WALL_KAPPA
unset CBS3D_WALL_B

# Native restart and exact global iteration continuation.
export CBS3D_RESTART=1
export CBS3D_RESTART_FORMAT=native
export CBS3D_RESTART_ROOT="$RESTART_ROOT"
export CBS3D_CONTINUATION_ITERATIONS=100000

# This continuation must leave native checkpoints behind.  Unlike the first
# 100k segment, we are not giving future us the opportunity to discover that a
# perfectly good final state was never checkpointed.
export CBS3D_CHECKPOINT_EVERY=10000
export CBS3D_CHECKPOINT_ROOT="$CHECKPOINT_ROOT"

export CBS3D_VERBOSE=1
unset CBS3D_ALL_RANK_OUTPUT

[[ -x "$EXE" ]] || { echo "FATAL: executable missing: $EXE" >&2; exit 1; }
[[ -f "$ANALYZER" ]] || { echo "FATAL: analyzer missing: $ANALYZER" >&2; exit 1; }
[[ -d "$SOURCE_PART" ]] || { echo "FATAL: source staged partition missing: $SOURCE_PART" >&2; exit 1; }
[[ -f "$RESTART_ROOT/restart_manifest.txt" ]] || { echo "FATAL: restart manifest missing" >&2; exit 1; }
[[ "$SLURM_NTASKS" -eq 40 ]] || { echo "FATAL: expected 40 MPI ranks" >&2; exit 1; }

rank_restart_count=$(find "$RESTART_ROOT" -maxdepth 1 -type f -name 'flatplate_step_00100000_rank_*.cbsrst' | wc -l)
[[ "$rank_restart_count" -eq 40 ]] || { echo "FATAL: expected 40 native restart rank files, found $rank_restart_count" >&2; exit 1; }

grep -q '^case_name flatplate$' "$RESTART_ROOT/restart_manifest.txt" || { echo "FATAL: restart case is not flatplate" >&2; exit 1; }
grep -q '^completed_iteration 100000$' "$RESTART_ROOT/restart_manifest.txt" || { echo "FATAL: restart is not iteration 100000" >&2; exit 1; }
grep -q '^mpi_size 40$' "$RESTART_ROOT/restart_manifest.txt" || { echo "FATAL: restart is not P40" >&2; exit 1; }

mkdir -p "$CASE/jobs/logs" "$JOB_ROOT" "$RUN_DIR"
rm -rf "$STAGED"
mkdir -p "$STAGED"
cp -a "$SOURCE_PART/." "$STAGED/"

rank_dirs=$(find "$STAGED" -maxdepth 1 -type d -name 'rank_[0-9][0-9][0-9][0-9]' | wc -l)
[[ "$rank_dirs" -eq 40 ]] || { echo "FATAL: expected 40 staged rank directories, found $rank_dirs" >&2; exit 1; }

printf '============================================================\n'
printf 'CBS3D SA TMR FLAT-PLATE NATIVE RESTART CONTINUATION\n'
printf '============================================================\n'
printf 'date                  = %s\n' "$(date -Is)"
printf 'host                  = %s\n' "$(hostname)"
printf 'job                   = %s\n' "$SLURM_JOB_ID"
printf 'source job            = %s\n' "$SOURCE_JOB_ID"
printf 'MPI ranks             = %s\n' "$SLURM_NTASKS"
printf 'restart root          = %s\n' "$RESTART_ROOT"
printf 'source partition      = %s\n' "$SOURCE_PART"
printf 'staged partition      = %s\n' "$STAGED"
printf 'run directory         = %s\n' "$RUN_DIR"
printf 'wall model            = OFF\n'
printf 'source iteration      = 100000\n'
printf 'additional iterations = 100000\n'
printf 'target global iter    = 200000\n'
printf 'checkpoint every      = 10000\n'
printf '============================================================\n'
cat "$RESTART_ROOT/restart_manifest.txt"
printf '============================================================\n'

cd "$RUN_DIR"

set +e
srun \
    --mpi=pmix \
    --unbuffered \
    --nodes=1 \
    --ntasks=40 \
    --ntasks-per-node=40 \
    --cpus-per-task=1 \
    --distribution=block:block \
    --cpu-bind=cores \
    --kill-on-bad-exit=1 \
    stdbuf -oL -eL \
    "$EXE" flatplate "$STAGED" \
    > solver.out \
    2> solver.err
RC=$?
set -e

echo
echo "===== FINAL TURBULENCE / FLOW EVIDENCE ====="
grep -E \
    'Restart loaded|Continuation segment|Iteration|RelU|RelV|RelW|RelP|RelSA|NuTilde|NuT=|MuT|MuEff|DivRMS|CG=|Run complete|stop reason|FATAL|ERROR' \
    solver.out | tail -n 280 || true

echo
echo "===== LAST 180 SOLVER LINES ====="
tail -n 180 solver.out || true

echo
echo "===== STDERR ====="
cat solver.err || true

echo
echo "RETURN CODE = $RC"

if [[ "$RC" -ne 0 ]]; then
    echo "CBS3D SA TMR RESTART CONTINUATION: FAILED" >&2
    exit "$RC"
fi

if grep -Eqi 'nan|FATAL|ERROR:' solver.out solver.err; then
    echo "FATAL: non-finite/fatal diagnostic detected in continuation logs." >&2
    exit 3
fi

if grep -q 'SA production wall treatment: ON' solver.out; then
    echo "FATAL: wall treatment was active in the wall-resolved benchmark." >&2
    exit 4
fi

FINAL_PVTU=$(find "$RUN_DIR/output/flatplate" -maxdepth 1 -type f -name '*.pvtu' 2>/dev/null | sort | tail -1 || true)
printf 'FINAL_PVTU=%s\n' "$FINAL_PVTU"
printf 'BENCHMARK_ROOT=%s\n' "$JOB_ROOT"

if [[ -z "$FINAL_PVTU" ]]; then
    echo "FATAL: continuation completed but no distributed VTU snapshot was found." >&2
    exit 5
fi

echo
echo "===== DISTRIBUTED TMR QUANTITATIVE DIAGNOSTICS ====="
python3 "$ANALYZER" \
    "$STAGED" \
    "$RUN_DIR/output/flatplate" \
    --station 0.97008 \
    --re-x1 5.0e6 \
    --u-inf 1.0 \
    --kappa 0.41 \
    | tee "$JOB_ROOT/tmr_diagnostics.json"

echo "TMR_DIAGNOSTICS=$JOB_ROOT/tmr_diagnostics.json"
echo "CBS3D SA TMR RESTART CONTINUATION: COMPLETED"
