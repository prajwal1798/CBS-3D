#!/usr/bin/env bash
# Run ONE matched-restart csafm sensitivity case for the wall-resolved NASA/TMR
# flat plate.  No solver equations are changed.
#
# Usage:
#   bash jobs/sunbird/run_flatplate_sa_csafm_single_restart.sh 0.50 [SOURCE_JOB_ID]
#   bash jobs/sunbird/run_flatplate_sa_csafm_single_restart.sh 0.25 [SOURCE_JOB_ID]
#   bash jobs/sunbird/run_flatplate_sa_csafm_single_restart.sh 0.10 [SOURCE_JOB_ID]
#
# Each branch starts from the same exact global-iteration-100000 P40 field.
# Continuation lengths keep N_cont*csafm = 5000:
#   0.50 -> 10000 iterations -> global 110000
#   0.25 -> 20000 iterations -> global 120000
#   0.10 -> 50000 iterations -> global 150000

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
RUNS=$CASE/runs/sa_tmr_wall_resolved
LOGS=$CASE/jobs/logs
GENERATED=$CASE/jobs/generated/csafm_single
EXE=$CODE/build/mpi-release/cbs3d_parallel
CONVERTER=$CODE/tools/flatplate_vtu_to_native_restart.py
ANALYZER=$CODE/examples/Flat_Plate_Turbulent/analyze_distributed_tmr.py

GCC_ROOT=/lustrehome/apps/compilers/gnu/12.1.0
MPI_ROOT=/apps/libraries/openmpi/4.1.6/el7/AVX512/gnu-12.1
export PETSC_DIR=/home/s.2337862/software/petsc-gcc12-install

CSAFM=${1:-}
SOURCE_JOB_ID=${2:-}

case "$CSAFM" in
    0.50|0.5)
        CSAFM=0.50
        TAG=050
        CONT_ITERS=10000
        TARGET_GLOBAL=110000
        ;;
    0.25)
        TAG=025
        CONT_ITERS=20000
        TARGET_GLOBAL=120000
        ;;
    0.10|0.1)
        CSAFM=0.10
        TAG=010
        CONT_ITERS=50000
        TARGET_GLOBAL=150000
        ;;
    *)
        echo "Usage: $0 {0.50|0.25|0.10} [SOURCE_JOB_ID]" >&2
        exit 2
        ;;
esac

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

if [[ -z "$SOURCE_JOB_ID" ]]; then
    SOURCE_JOB_ID=$(
        find "$RUNS" -maxdepth 1 -type d -name 'job_*' -printf '%f\n' 2>/dev/null \
        | sed 's/^job_//' \
        | sort -nr \
        | while read -r job; do
            out="$RUNS/job_${job}/run/solver.out"
            if [[ -f "$out" ]] \
               && grep -q '^Run complete$' "$out" \
               && grep -q '^Iteration 100000/100000 ' "$out"; then
                echo "$job"
                break
            fi
          done
    )
fi

[[ -n "$SOURCE_JOB_ID" ]] || {
    echo "FATAL: no completed exact-100k flat-plate source job found" >&2
    exit 1
}

SOURCE_JOB=$RUNS/job_${SOURCE_JOB_ID}
SOURCE_RUN=$SOURCE_JOB/run
SOURCE_OUT=$SOURCE_RUN/output/flatplate
SOURCE_PART=$SOURCE_JOB/partitions
SOURCE_SOLVER_OUT=$SOURCE_RUN/solver.out
SOURCE_STEP=100000
SOURCE_STEP_TAG=step_00100000
RESTART_DIR=$SOURCE_OUT/restart_from_vtu/csafm_single_${SOURCE_STEP_TAG}

[[ -d "$SOURCE_PART" ]] || { echo "FATAL: source partition missing: $SOURCE_PART" >&2; exit 1; }
[[ -d "$SOURCE_OUT" ]] || { echo "FATAL: source output missing: $SOURCE_OUT" >&2; exit 1; }
[[ -f "$SOURCE_SOLVER_OUT" ]] || { echo "FATAL: source solver.out missing" >&2; exit 1; }
[[ -f "$CONVERTER" ]] || { echo "FATAL: restart converter missing: $CONVERTER" >&2; exit 1; }
[[ -f "$ANALYZER" ]] || { echo "FATAL: TMR analyzer missing: $ANALYZER" >&2; exit 1; }

grep -q '^Run complete$' "$SOURCE_SOLVER_OUT" || {
    echo "FATAL: source run did not complete" >&2
    exit 1
}
grep -q '^Iteration 100000/100000 ' "$SOURCE_SOLVER_OUT" || {
    echo "FATAL: source run is not exact global iteration 100000" >&2
    exit 1
}

piece_count=$(find "$SOURCE_OUT" -maxdepth 1 -type f -name "flatplate_${SOURCE_STEP_TAG}_rank_*.vtu" | wc -l)
[[ "$piece_count" -eq 40 ]] || {
    echo "FATAL: expected 40 source VTU pieces at 100000; found $piece_count" >&2
    exit 1
}
[[ -f "$SOURCE_OUT/flatplate_${SOURCE_STEP_TAG}.pvtu" ]] || {
    echo "FATAL: source final PVTU missing" >&2
    exit 1
}

DT=$(python3 - "$SOURCE_SOLVER_OUT" <<'PY'
from __future__ import print_function
import re
import sys

pat = re.compile(r"^Iteration\s+100000/100000\s+dt=([^\s]+)")
value = None
with open(sys.argv[1], "r") as handle:
    for line in handle:
        m = pat.search(line)
        if m:
            value = m.group(1)
if value is None:
    raise SystemExit("FATAL: final source dt not found")
print(value)
PY
)

printf 'source job       = %s\n' "$SOURCE_JOB_ID"
printf 'source final dt  = %s\n' "$DT"
printf 'requested csafm  = %s\n' "$CSAFM"
printf 'continuation     = %s\n' "$CONT_ITERS"
printf 'target global    = %s\n' "$TARGET_GLOBAL"

printf '\n===== CREATE / VERIFY COMMON NATIVE RESTART =====\n'
if [[ ! -f "$RESTART_DIR/restart_manifest.txt" ]]; then
    rm -rf "$RESTART_DIR"
    python3 "$CONVERTER" \
        --input-dir "$SOURCE_OUT" \
        --output-dir "$RESTART_DIR" \
        --case flatplate \
        --iteration "$SOURCE_STEP" \
        --mpi-size 40 \
        --physical-time 0.0 \
        --time-step "$DT" \
        --temperature-enabled 0 \
        --turbulence-enabled 1
fi

[[ -f "$RESTART_DIR/restart_manifest.txt" ]] || {
    echo "FATAL: restart manifest missing" >&2
    exit 1
}
restart_count=$(find "$RESTART_DIR" -maxdepth 1 -type f -name 'flatplate_step_00100000_rank_*.cbsrst' | wc -l)
[[ "$restart_count" -eq 40 ]] || {
    echo "FATAL: expected 40 native restart files; found $restart_count" >&2
    exit 1
}
grep -q '^completed_iteration 100000$' "$RESTART_DIR/restart_manifest.txt" || {
    echo "FATAL: restart manifest is not iteration 100000" >&2
    exit 1
}
grep -q '^mpi_size 40$' "$RESTART_DIR/restart_manifest.txt" || {
    echo "FATAL: restart manifest is not P40" >&2
    exit 1
}

printf '\n===== BUILD MPI/PETSC RELEASE =====\n'
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

mkdir -p "$LOGS" "$GENERATED"
JOB_SCRIPT=$GENERATED/flatplate_sa_csafm_single_${HEAD}.slurm

cat > "$JOB_SCRIPT" <<'SLURM'
#!/bin/bash -l
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/s.2337862/CBS3D/cases/flatplate/jobs/logs/%x_%j.out
#SBATCH --error=/scratch/s.2337862/CBS3D/cases/flatplate/jobs/logs/%x_%j.err

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
EXE=$CODE/build/mpi-release/cbs3d_parallel
ANALYZER=$CODE/examples/Flat_Plate_Turbulent/analyze_distributed_tmr.py

: "${SOURCE_JOB_ID:?missing SOURCE_JOB_ID}"
: "${RESTART_ROOT:?missing RESTART_ROOT}"
: "${CSAFM:?missing CSAFM}"
: "${CONT_ITERS:?missing CONT_ITERS}"
: "${TARGET_GLOBAL:?missing TARGET_GLOBAL}"
: "${TAG:?missing TAG}"

SOURCE_JOB=$CASE/runs/sa_tmr_wall_resolved/job_${SOURCE_JOB_ID}
SOURCE_PART=$SOURCE_JOB/partitions
RUN_ROOT=$CASE/runs/sa_tmr_csafm_single
JOB_ROOT=$RUN_ROOT/csafm_${TAG}_job_${SLURM_JOB_ID}
STAGED=$JOB_ROOT/partitions
RUN_DIR=$JOB_ROOT/run

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
export OMPI_MCA_btl="^openib"

unset CBS3D_SA_WALL_TREATMENT
unset CBS3D_CHT_WALL_TREATMENT
unset CBS3D_WALL_KAPPA
unset CBS3D_WALL_B

export CBS3D_RESTART=1
export CBS3D_RESTART_FORMAT=native
export CBS3D_RESTART_ROOT="$RESTART_ROOT"
export CBS3D_CONTINUATION_ITERATIONS="$CONT_ITERS"
export CBS3D_VERBOSE=1
unset CBS3D_ALL_RANK_OUTPUT
unset CBS3D_CHECKPOINT_EVERY
unset CBS3D_CHECKPOINT_ROOT

[[ -x "$EXE" ]] || { echo "FATAL: executable missing: $EXE" >&2; exit 1; }
[[ -f "$ANALYZER" ]] || { echo "FATAL: analyzer missing: $ANALYZER" >&2; exit 1; }
[[ -d "$SOURCE_PART" ]] || { echo "FATAL: source partition missing: $SOURCE_PART" >&2; exit 1; }
[[ -f "$RESTART_ROOT/restart_manifest.txt" ]] || { echo "FATAL: restart manifest missing" >&2; exit 1; }
[[ "$SLURM_NTASKS" -eq 40 ]] || { echo "FATAL: expected 40 MPI ranks" >&2; exit 1; }

grep -q '^completed_iteration 100000$' "$RESTART_ROOT/restart_manifest.txt" || {
    echo "FATAL: restart source is not global iteration 100000" >&2
    exit 1
}
grep -q '^mpi_size 40$' "$RESTART_ROOT/restart_manifest.txt" || {
    echo "FATAL: restart source is not P40" >&2
    exit 1
}

mkdir -p "$JOB_ROOT" "$RUN_DIR"
rm -rf "$STAGED"
mkdir -p "$STAGED"
cp -a "$SOURCE_PART/." "$STAGED/"

# Patch csafm by the same positional contract used by MeshIO::readParameterFile.
# The solver does not inspect the text of label lines.  It consumes alternating
# label/data records.  CSAFM/THETA is the sixth fixed data block, i.e. active
# (nonblank, noncomment) line index 11 in zero-based Python indexing.
python3 - "$STAGED" "$CSAFM" <<'PY'
from __future__ import print_function
import math
import os
import sys

root = sys.argv[1]
new_csafm = float(sys.argv[2])
if new_csafm not in (0.5, 0.25, 0.1):
    raise SystemExit("FATAL: unexpected csafm {}".format(new_csafm))

rank_dirs = sorted(
    os.path.join(root, name)
    for name in os.listdir(root)
    if name.startswith("rank_") and os.path.isdir(os.path.join(root, name))
)
if len(rank_dirs) != 40:
    raise SystemExit("FATAL: expected 40 rank directories; found {}".format(len(rank_dirs)))

old_values = []
theta_values = []

for rank_dir in rank_dirs:
    pars = sorted(
        os.path.join(rank_dir, name)
        for name in os.listdir(rank_dir)
        if name.endswith(".par")
    )
    if len(pars) != 1:
        raise SystemExit("FATAL: expected one .par in {}".format(rank_dir))
    path = pars[0]
    with open(path, "r") as handle:
        lines = handle.read().splitlines()

    active = [
        i for i, line in enumerate(lines)
        if line.strip() and not line.lstrip().startswith(("#", "!"))
    ]
    if len(active) < 14:
        raise SystemExit("FATAL: too few active .par lines in {}".format(path))

    # Structural guards around the fixed early-file layout read by MeshIO.
    time_tokens = lines[active[9]].split()
    cs_tokens = lines[active[11]].split()
    cbs_tokens = lines[active[13]].split()
    if len(time_tokens) < 5:
        raise SystemExit("FATAL: timestep row at active[9] malformed in {}".format(path))
    if len(cs_tokens) != 2:
        raise SystemExit("FATAL: CSAFM/THETA row at active[11] malformed in {}: {!r}".format(path, lines[active[11]]))
    if len(cbs_tokens) < 10:
        raise SystemExit("FATAL: CBS timestep row at active[13] malformed in {}".format(path))

    try:
        old_csafm = float(cs_tokens[0])
        theta = float(cs_tokens[1])
        [float(v) for v in time_tokens[:5]]
        [float(v) for v in cbs_tokens[:10]]
    except ValueError:
        raise SystemExit("FATAL: numeric parse failed around CSAFM block in {}".format(path))

    if not (old_csafm > 0.0 and math.isfinite(old_csafm)):
        raise SystemExit("FATAL: invalid source csafm in {}".format(path))
    if not (0.5 <= theta <= 1.0 and math.isfinite(theta)):
        raise SystemExit("FATAL: invalid theta in {}".format(path))

    lines[active[11]] = "{:.17g} {:.17g}".format(new_csafm, theta)
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")

    old_values.append(old_csafm)
    theta_values.append(theta)

old_unique = sorted(set(round(v, 15) for v in old_values))
theta_unique = sorted(set(round(v, 15) for v in theta_values))
if len(old_unique) != 1:
    raise SystemExit("FATAL: source partitions disagree on csafm: {}".format(old_unique))
if len(theta_unique) != 1:
    raise SystemExit("FATAL: source partitions disagree on theta: {}".format(theta_unique))

print("source csafm      = {}".format(old_unique[0]))
print("requested csafm   = {:.17g}".format(new_csafm))
print("theta             = {}".format(theta_unique[0]))
print("patched .par files= 40")
PY

printf '============================================================\n'
printf 'CBS3D SINGLE CSAFM MATCHED-RESTART CASE\n'
printf '============================================================\n'
printf 'job            = %s\n' "$SLURM_JOB_ID"
printf 'source job     = %s\n' "$SOURCE_JOB_ID"
printf 'csafm          = %s\n' "$CSAFM"
printf 'restart global = 100000\n'
printf 'continuation   = %s\n' "$CONT_ITERS"
printf 'target global  = %s\n' "$TARGET_GLOBAL"
printf 'wall model     = OFF\n'
printf 'run dir        = %s\n' "$RUN_DIR"
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

printf '\n===== FINAL ITERATION =====\n'
grep '^Iteration ' solver.out | tail -n 1 || true
printf '\n===== STDERR =====\n'
cat solver.err || true
printf '\nRETURN CODE = %s\n' "$RC"

if [[ "$RC" -ne 0 ]]; then
    echo "SINGLE CSAFM CASE FAILED" >&2
    exit "$RC"
fi
if grep -Eqi 'nan|FATAL|ERROR:' solver.out solver.err; then
    echo "FATAL: non-finite/fatal diagnostic detected" >&2
    exit 3
fi
if grep -q 'SA production wall treatment: ON' solver.out; then
    echo "FATAL: wall treatment active in wall-resolved case" >&2
    exit 4
fi
if ! grep -q '^Run complete$' solver.out; then
    echo "FATAL: solver did not report Run complete" >&2
    exit 5
fi
if ! grep -q "^Iteration ${TARGET_GLOBAL}/${TARGET_GLOBAL} " solver.out; then
    echo "FATAL: requested target global iteration ${TARGET_GLOBAL} was not reached" >&2
    exit 6
fi

FINAL_PVTU=$(find "$RUN_DIR/output/flatplate" -maxdepth 1 -type f -name '*.pvtu' 2>/dev/null | sort | tail -1 || true)
[[ -n "$FINAL_PVTU" ]] || { echo "FATAL: final PVTU missing" >&2; exit 7; }
printf 'FINAL_PVTU=%s\n' "$FINAL_PVTU"

printf '\n===== TMR QUANTITATIVE DIAGNOSTICS =====\n'
python3 "$ANALYZER" \
    "$STAGED" \
    "$RUN_DIR/output/flatplate" \
    --station 0.97008 \
    --re-x1 5.0e6 \
    --u-inf 1.0 \
    --kappa 0.41 \
    | tee "$JOB_ROOT/tmr_diagnostics.json"

printf '\nCSAFM_RESULT_ROOT=%s\n' "$JOB_ROOT"
printf 'CSAFM_DIAGNOSTICS=%s/tmr_diagnostics.json\n' "$JOB_ROOT"
printf 'SINGLE CSAFM CASE: COMPLETED\n'
SLURM

chmod 700 "$JOB_SCRIPT"

JOB_NAME=fp_c${TAG}_one
printf '\n===== SUBMIT ONE CASE ONLY =====\n'
JOB_ID=$(sbatch \
    --parsable \
    --job-name="$JOB_NAME" \
    --export=ALL,SOURCE_JOB_ID="$SOURCE_JOB_ID",RESTART_ROOT="$RESTART_DIR",CSAFM="$CSAFM",CONT_ITERS="$CONT_ITERS",TARGET_GLOBAL="$TARGET_GLOBAL",TAG="$TAG" \
    "$JOB_SCRIPT")

printf 'JOBID         = %s\n' "$JOB_ID"
printf 'JOB NAME      = %s\n' "$JOB_NAME"
printf 'CSAFM         = %s\n' "$CSAFM"
printf 'SOURCE JOB    = %s\n' "$SOURCE_JOB_ID"
printf 'TARGET GLOBAL = %s\n' "$TARGET_GLOBAL"
printf 'OUT           = %s/%s_%s.out\n' "$LOGS" "$JOB_NAME" "$JOB_ID"
printf 'ERR           = %s/%s_%s.err\n' "$LOGS" "$JOB_NAME" "$JOB_ID"
printf 'RUN ROOT      = %s/runs/sa_tmr_csafm_single/csafm_%s_job_%s\n' "$CASE" "$TAG" "$JOB_ID"
printf '\nQueue:\n'
squeue -j "$JOB_ID" || true
printf '\nLive solver output after startup:\n'
printf 'tail -F %s/runs/sa_tmr_csafm_single/csafm_%s_job_%s/run/solver.out\n' "$CASE" "$TAG" "$JOB_ID"
