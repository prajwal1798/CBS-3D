#!/usr/bin/env bash
# Matched-restart csafm sensitivity sweep for the wall-resolved NASA/TMR
# flat-plate SA benchmark.  No solver equations are changed by this script.
#
# All branches start from the exact same completed 100000-iteration P40 field.
# Continuation lengths are chosen so that N_cont * csafm = 5000 for each case:
#
#     csafm=0.50 -> 10000 iterations
#     csafm=0.25 -> 20000 iterations
#     csafm=0.10 -> 50000 iterations
#
# This is a nominal equal-pseudo-time comparison.  LTS evolves with the state,
# so it is not an exact integral of every local dt, but it avoids the much worse
# fixed-iteration comparison where smaller csafm simply advances less pseudo-time.

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
RUNS=$CASE/runs/sa_tmr_wall_resolved
LOGS=$CASE/jobs/logs
GENERATED=$CASE/jobs/generated/csafm_sweep
EXE=$CODE/build/mpi-release/cbs3d_parallel
CONVERTER=$CODE/tools/flatplate_vtu_to_native_restart.py
ANALYZER=$CODE/examples/Flat_Plate_Turbulent/analyze_distributed_tmr.py

GCC_ROOT=/lustrehome/apps/compilers/gnu/12.1.0
MPI_ROOT=/apps/libraries/openmpi/4.1.6/el7/AVX512/gnu-12.1
export PETSC_DIR=/home/s.2337862/software/petsc-gcc12-install

SOURCE_JOB_ID=${1:-}

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
TAG=step_00100000
RESTART_DIR=$SOURCE_OUT/restart_from_vtu/csafm_sweep_${TAG}

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
    echo "FATAL: source run is not an exact 100000-iteration run" >&2
    exit 1
}

piece_count=$(find "$SOURCE_OUT" -maxdepth 1 -type f -name "flatplate_${TAG}_rank_*.vtu" | wc -l)
[[ "$piece_count" -eq 40 ]] || {
    echo "FATAL: expected 40 source VTU pieces at iteration 100000; found $piece_count" >&2
    exit 1
}
[[ -f "$SOURCE_OUT/flatplate_${TAG}.pvtu" ]] || {
    echo "FATAL: source final PVTU missing" >&2
    exit 1
}

DT=$(python3 - "$SOURCE_SOLVER_OUT" <<'PY'
from __future__ import print_function
import re
import sys

pattern = re.compile(r"^Iteration\s+100000/100000\s+dt=([^\s]+)")
value = None
with open(sys.argv[1], "r") as handle:
    for line in handle:
        match = pattern.search(line)
        if match:
            value = match.group(1)
if value is None:
    raise SystemExit("FATAL: final source dt not found")
print(value)
PY
)

printf 'source job       = %s\n' "$SOURCE_JOB_ID"
printf 'source final dt  = %s\n' "$DT"
printf 'source VTU pieces= %s\n' "$piece_count"

printf '\n===== CREATE COMMON NATIVE RESTART =====\n'
rm -rf "$RESTART_DIR"
python3 "$CONVERTER" \
    --input-dir "$SOURCE_OUT" \
    --output-dir "$RESTART_DIR" \
    --case flatplate \
    --iteration 100000 \
    --mpi-size 40 \
    --physical-time 0.0 \
    --time-step "$DT" \
    --temperature-enabled 0 \
    --turbulence-enabled 1

[[ -f "$RESTART_DIR/restart_manifest.txt" ]] || {
    echo "FATAL: restart manifest missing after conversion" >&2
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
JOB_SCRIPT=$GENERATED/flatplate_sa_csafm_restart_sweep_${HEAD}.slurm

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
: "${SWEEP_TAG:?missing SWEEP_TAG}"

SOURCE_JOB=$CASE/runs/sa_tmr_wall_resolved/job_${SOURCE_JOB_ID}
SOURCE_PART=$SOURCE_JOB/partitions
RUN_ROOT=$CASE/runs/sa_tmr_csafm_sweep
JOB_ROOT=$RUN_ROOT/${SWEEP_TAG}_job_${SLURM_JOB_ID}
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

python3 - "$STAGED" "$CSAFM" "$TARGET_GLOBAL" <<'PY'
from __future__ import print_function
import os
import sys

root = sys.argv[1]
csafm = float(sys.argv[2])
target_global = int(sys.argv[3])

if csafm not in (0.5, 0.25, 0.1):
    raise SystemExit("FATAL: unexpected csafm {}".format(csafm))
if target_global <= 100000:
    raise SystemExit("FATAL: invalid target global iteration")

rank_dirs = sorted(
    os.path.join(root, name)
    for name in os.listdir(root)
    if name.startswith("rank_") and os.path.isdir(os.path.join(root, name))
)
if len(rank_dirs) != 40:
    raise SystemExit("FATAL: expected 40 staged rank directories; found {}".format(len(rank_dirs)))


def active_indices(lines):
    return [
        i for i, line in enumerate(lines)
        if line.strip() and not line.lstrip().startswith(("#", "!"))
    ]


def unique_label_data(lines, predicate, description):
    active = active_indices(lines)
    hits = [p for p, i in enumerate(active) if predicate(lines[i].strip().lower())]
    if len(hits) != 1:
        raise SystemExit("FATAL: expected one {} label; found {}".format(description, len(hits)))
    pos = hits[0]
    if pos + 1 >= len(active):
        raise SystemExit("FATAL: missing data after {} label".format(description))
    return active[pos + 1]

seen = []
for rank_dir in rank_dirs:
    par_files = sorted(
        os.path.join(rank_dir, name)
        for name in os.listdir(rank_dir)
        if name.endswith(".par")
    )
    if len(par_files) != 1:
        raise SystemExit("FATAL: expected exactly one .par in {}".format(rank_dir))
    path = par_files[0]
    with open(path, "r") as handle:
        lines = handle.read().splitlines()

    cs_i = unique_label_data(
        lines,
        lambda s: s == "csafm theta",
        "csafm theta",
    )
    vals = lines[cs_i].split()
    if len(vals) < 2:
        raise SystemExit("FATAL: malformed csafm/theta row in {}".format(path))
    old_csafm = float(vals[0])
    theta = vals[1]
    vals[0] = "{:.17g}".format(csafm)
    lines[cs_i] = "  ".join(vals)

    out_i = unique_label_data(
        lines,
        lambda s: s.startswith("cbs3d output/monitor controls:"),
        "output/monitor",
    )
    out = lines[out_i].split()
    if len(out) < 9:
        raise SystemExit("FATAL: malformed output/monitor row in {}".format(path))
    out[0] = "1"       # residual logging on
    out[1] = "10"      # residual every 10
    out[2] = "1000"    # console every 1000
    out[3] = "0"       # no live plotting
    out[4] = "1"       # VTU enabled
    out[5] = "10000"   # VTU every 10000 global iterations
    out[8] = str(target_global)  # forbid convergence stop before target
    lines[out_i] = "  ".join(out)

    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    seen.append((old_csafm, theta))

old_values = sorted(set(round(v[0], 15) for v in seen))
theta_values = sorted(set(v[1] for v in seen))
if len(old_values) != 1:
    raise SystemExit("FATAL: source partitions disagree on csafm: {}".format(old_values))
if len(theta_values) != 1:
    raise SystemExit("FATAL: source partitions disagree on theta: {}".format(theta_values))

print("source csafm = {}".format(old_values[0]))
print("sweep csafm  = {:.17g}".format(csafm))
print("theta        = {}".format(theta_values[0]))
print("target iter  = {}".format(target_global))
print("patched rank parameter files = 40")
PY

printf '============================================================\n'
printf 'CBS3D MATCHED-RESTART CSAFM SWEEP\n'
printf '============================================================\n'
printf 'job              = %s\n' "$SLURM_JOB_ID"
printf 'source job       = %s\n' "$SOURCE_JOB_ID"
printf 'csafm            = %s\n' "$CSAFM"
printf 'restart global   = 100000\n'
printf 'continuation     = %s\n' "$CONT_ITERS"
printf 'target global    = %s\n' "$TARGET_GLOBAL"
printf 'nominal N*csafm = 5000\n'
printf 'wall model       = OFF\n'
printf 'run dir          = %s\n' "$RUN_DIR"
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
    echo "CSAFM SWEEP BRANCH FAILED" >&2
    exit "$RC"
fi
if grep -Eqi 'nan|FATAL|ERROR:' solver.out solver.err; then
    echo "FATAL: non-finite/fatal diagnostic detected" >&2
    exit 3
fi
if grep -q 'SA production wall treatment: ON' solver.out; then
    echo "FATAL: wall treatment active in wall-resolved sweep" >&2
    exit 4
fi

FINAL_PVTU=$(find "$RUN_DIR/output/flatplate" -maxdepth 1 -type f -name '*.pvtu' | sort | tail -1 || true)
[[ -n "$FINAL_PVTU" ]] || { echo "FATAL: final PVTU missing" >&2; exit 5; }
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
SLURM

chmod 700 "$JOB_SCRIPT"

printf '\n===== SUBMIT MATCHED-RESTART SWEEP =====\n'
printf 'Common source job    : %s\n' "$SOURCE_JOB_ID"
printf 'Common restart root  : %s\n' "$RESTART_DIR"
printf 'Generated job script : %s\n' "$JOB_SCRIPT"

submit_case() {
    local csafm=$1
    local tag=$2
    local cont=$3
    local target=$4
    local name="fp_c${tag}"
    local jid
    jid=$(sbatch \
        --parsable \
        --job-name="$name" \
        --export=ALL,SOURCE_JOB_ID="$SOURCE_JOB_ID",RESTART_ROOT="$RESTART_DIR",CSAFM="$csafm",CONT_ITERS="$cont",TARGET_GLOBAL="$target",SWEEP_TAG="csafm_${tag}" \
        "$JOB_SCRIPT")
    printf 'csafm=%-5s  continuation=%-6s  target=%-6s  JOBID=%s  log=%s/%s_%s.out\n' \
        "$csafm" "$cont" "$target" "$jid" "$LOGS" "$name" "$jid"
}

submit_case 0.50 050 10000 110000
submit_case 0.25 025 20000 120000
submit_case 0.10 010 50000 150000

printf '\nQueue:\n'
squeue -u "$USER" -n fp_c050,fp_c025,fp_c010 || true

printf '\nAfter all three finish, compare diagnostics with:\n'
printf 'for f in %s/runs/sa_tmr_csafm_sweep/csafm_*/tmr_diagnostics.json; do echo "===== $f ====="; cat "$f"; done\n' "$CASE"
