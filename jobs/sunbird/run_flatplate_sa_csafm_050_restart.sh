#!/usr/bin/env bash
# Dedicated one-at-a-time csafm=0.50 continuation for the wall-resolved
# NASA/TMR flat-plate benchmark.  Starts from the exact P40 global-100000 field
# and advances 10000 segment-local iterations to global iteration 110000.

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
RUNS=$CASE/runs/sa_tmr_wall_resolved
LOGS=$CASE/jobs/logs
GENERATED=$CASE/jobs/generated/csafm_050
EXE=$CODE/build/mpi-release/cbs3d_parallel
CONVERTER=$CODE/tools/flatplate_vtu_to_native_restart.py
ANALYZER=$CODE/examples/Flat_Plate_Turbulent/analyze_distributed_tmr.py

GCC_ROOT=/lustrehome/apps/compilers/gnu/12.1.0
MPI_ROOT=/apps/libraries/openmpi/4.1.6/el7/AVX512/gnu-12.1
export PETSC_DIR=/home/s.2337862/software/petsc-gcc12-install

CSAFM=0.50
TAG=050
CONT_ITERS=10000
TARGET_GLOBAL=110000
TARGET_STEP_TAG=step_00110000
SOURCE_JOB_ID=${1:-}

cd "$CODE"

printf '===== SOURCE GATE =====\n'
BRANCH=$(git rev-parse --abbrev-ref HEAD)
HEAD=$(git rev-parse HEAD)
printf 'branch = %s\n' "$BRANCH"
printf 'HEAD   = %s\n' "$HEAD"
[[ "$BRANCH" == "integration/sa-mpi" ]] || { echo "FATAL: expected integration/sa-mpi" >&2; exit 1; }
[[ -z "$(git status --porcelain)" ]] || {
    echo "FATAL: Sunbird worktree is dirty" >&2
    git status --short
    exit 1
}

if [[ -z "$SOURCE_JOB_ID" ]]; then
    SOURCE_JOB_ID=$(
        find "$RUNS" -maxdepth 1 -type d -name 'job_*' -printf '%f\n' 2>/dev/null \
        | sed 's/^job_//' | sort -nr \
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
[[ -n "$SOURCE_JOB_ID" ]] || { echo "FATAL: no completed exact-100k source job found" >&2; exit 1; }

SOURCE_JOB=$RUNS/job_${SOURCE_JOB_ID}
SOURCE_RUN=$SOURCE_JOB/run
SOURCE_OUT=$SOURCE_RUN/output/flatplate
SOURCE_PART=$SOURCE_JOB/partitions
SOURCE_SOLVER_OUT=$SOURCE_RUN/solver.out
SOURCE_STEP_TAG=step_00100000
RESTART_DIR=$SOURCE_OUT/restart_from_vtu/csafm_050_${SOURCE_STEP_TAG}

[[ -d "$SOURCE_PART" ]] || { echo "FATAL: source partition missing: $SOURCE_PART" >&2; exit 1; }
[[ -d "$SOURCE_OUT" ]] || { echo "FATAL: source output missing: $SOURCE_OUT" >&2; exit 1; }
[[ -f "$SOURCE_SOLVER_OUT" ]] || { echo "FATAL: source solver.out missing" >&2; exit 1; }
[[ -f "$CONVERTER" ]] || { echo "FATAL: restart converter missing: $CONVERTER" >&2; exit 1; }
[[ -f "$ANALYZER" ]] || { echo "FATAL: analyzer missing: $ANALYZER" >&2; exit 1; }

grep -q '^Run complete$' "$SOURCE_SOLVER_OUT" || { echo "FATAL: source run did not complete" >&2; exit 1; }
grep -q '^Iteration 100000/100000 ' "$SOURCE_SOLVER_OUT" || { echo "FATAL: source is not exact 100k" >&2; exit 1; }

piece_count=$(find "$SOURCE_OUT" -maxdepth 1 -type f -name "flatplate_${SOURCE_STEP_TAG}_rank_*.vtu" | wc -l)
[[ "$piece_count" -eq 40 ]] || { echo "FATAL: expected 40 source VTUs; found $piece_count" >&2; exit 1; }
[[ -f "$SOURCE_OUT/flatplate_${SOURCE_STEP_TAG}.pvtu" ]] || { echo "FATAL: source PVTU missing" >&2; exit 1; }

DT=$(python3 - "$SOURCE_SOLVER_OUT" <<'PY'
from __future__ import print_function
import re, sys
pat = re.compile(r"^Iteration\s+100000/100000\s+dt=([^\s]+)")
value = None
with open(sys.argv[1], "r") as f:
    for line in f:
        m = pat.search(line)
        if m:
            value = m.group(1)
if value is None:
    raise SystemExit("FATAL: final source dt not found")
print(value)
PY
)

printf 'source job      = %s\n' "$SOURCE_JOB_ID"
printf 'source final dt = %s\n' "$DT"
printf 'csafm           = %s\n' "$CSAFM"
printf 'continuation    = %s local iterations\n' "$CONT_ITERS"
printf 'target global   = %s\n' "$TARGET_GLOBAL"

printf '\n===== CREATE / VERIFY NATIVE RESTART =====\n'
if [[ ! -f "$RESTART_DIR/restart_manifest.txt" ]]; then
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
fi
[[ -f "$RESTART_DIR/restart_manifest.txt" ]] || { echo "FATAL: restart manifest missing" >&2; exit 1; }
restart_count=$(find "$RESTART_DIR" -maxdepth 1 -type f -name 'flatplate_step_00100000_rank_*.cbsrst' | wc -l)
[[ "$restart_count" -eq 40 ]] || { echo "FATAL: expected 40 restart files; found $restart_count" >&2; exit 1; }
grep -q '^completed_iteration 100000$' "$RESTART_DIR/restart_manifest.txt" || { echo "FATAL: restart not global 100000" >&2; exit 1; }
grep -q '^mpi_size 40$' "$RESTART_DIR/restart_manifest.txt" || { echo "FATAL: restart not P40" >&2; exit 1; }

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
[[ -x "$EXE" ]] || { echo "FATAL: executable missing after build" >&2; exit 1; }

mkdir -p "$LOGS" "$GENERATED"
JOB_SCRIPT=$GENERATED/flatplate_sa_csafm_050_${HEAD}.slurm

cat > "$JOB_SCRIPT" <<'SLURM'
#!/bin/bash -l
#SBATCH --job-name=fp_c050_one
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/s.2337862/CBS3D/cases/flatplate/jobs/logs/fp_c050_one_%j.out
#SBATCH --error=/scratch/s.2337862/CBS3D/cases/flatplate/jobs/logs/fp_c050_one_%j.err

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
EXE=$CODE/build/mpi-release/cbs3d_parallel
ANALYZER=$CODE/examples/Flat_Plate_Turbulent/analyze_distributed_tmr.py

: "${SOURCE_JOB_ID:?missing SOURCE_JOB_ID}"
: "${RESTART_ROOT:?missing RESTART_ROOT}"

SOURCE_JOB=$CASE/runs/sa_tmr_wall_resolved/job_${SOURCE_JOB_ID}
SOURCE_PART=$SOURCE_JOB/partitions
RUN_ROOT=$CASE/runs/sa_tmr_csafm_single
JOB_ROOT=$RUN_ROOT/csafm_050_job_${SLURM_JOB_ID}
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
export CBS3D_CONTINUATION_ITERATIONS=10000
export CBS3D_VERBOSE=1
unset CBS3D_ALL_RANK_OUTPUT
unset CBS3D_CHECKPOINT_EVERY
unset CBS3D_CHECKPOINT_ROOT

[[ -x "$EXE" ]] || { echo "FATAL: executable missing" >&2; exit 1; }
[[ -f "$ANALYZER" ]] || { echo "FATAL: analyzer missing" >&2; exit 1; }
[[ -d "$SOURCE_PART" ]] || { echo "FATAL: source partitions missing" >&2; exit 1; }
[[ -f "$RESTART_ROOT/restart_manifest.txt" ]] || { echo "FATAL: restart manifest missing" >&2; exit 1; }
[[ "$SLURM_NTASKS" -eq 40 ]] || { echo "FATAL: expected 40 MPI ranks" >&2; exit 1; }

mkdir -p "$JOB_ROOT" "$RUN_DIR"
rm -rf "$STAGED"
mkdir -p "$STAGED"
cp -a "$SOURCE_PART/." "$STAGED/"

# Mirror MeshIO::readParameterFile's fixed early-file positional contract.
python3 - "$STAGED" <<'PY'
from __future__ import print_function
import math, os, sys
root = sys.argv[1]
rank_dirs = sorted(os.path.join(root, n) for n in os.listdir(root)
                   if n.startswith("rank_") and os.path.isdir(os.path.join(root, n)))
if len(rank_dirs) != 40:
    raise SystemExit("FATAL: expected 40 staged ranks; found {}".format(len(rank_dirs)))
old_values = []
theta_values = []
for rank_dir in rank_dirs:
    pars = sorted(os.path.join(rank_dir, n) for n in os.listdir(rank_dir) if n.endswith(".par"))
    if len(pars) != 1:
        raise SystemExit("FATAL: expected one .par in {}".format(rank_dir))
    path = pars[0]
    with open(path, "r") as f:
        lines = f.read().splitlines()
    active = [i for i, line in enumerate(lines)
              if line.strip() and not line.lstrip().startswith(("#", "!"))]
    if len(active) < 14:
        raise SystemExit("FATAL: too few active .par lines in {}".format(path))
    time_tokens = lines[active[9]].split()
    cs_tokens = lines[active[11]].split()
    cbs_tokens = lines[active[13]].split()
    if len(time_tokens) < 5 or len(cs_tokens) != 2 or len(cbs_tokens) < 10:
        raise SystemExit("FATAL: early .par positional structure mismatch in {}".format(path))
    try:
        [float(v) for v in time_tokens[:5]]
        old_csafm = float(cs_tokens[0])
        theta = float(cs_tokens[1])
        [float(v) for v in cbs_tokens[:10]]
    except ValueError:
        raise SystemExit("FATAL: numeric parse failed around CSAFM block in {}".format(path))
    if not (old_csafm > 0.0 and math.isfinite(old_csafm)):
        raise SystemExit("FATAL: invalid source CSAFM in {}".format(path))
    if not (0.5 <= theta <= 1.0 and math.isfinite(theta)):
        raise SystemExit("FATAL: invalid THETA in {}".format(path))
    lines[active[11]] = "0.5 {:.17g}".format(theta)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    old_values.append(old_csafm)
    theta_values.append(theta)
if len(set(round(v, 15) for v in old_values)) != 1:
    raise SystemExit("FATAL: source ranks disagree on CSAFM")
if len(set(round(v, 15) for v in theta_values)) != 1:
    raise SystemExit("FATAL: source ranks disagree on THETA")
print("source CSAFM      = {}".format(old_values[0]))
print("requested CSAFM   = 0.5")
print("theta             = {}".format(theta_values[0]))
print("patched .par files= 40")
PY

printf '============================================================\n'
printf 'CBS3D CSAFM 0.50 SINGLE RESTART\n'
printf '============================================================\n'
printf 'job            = %s\n' "$SLURM_JOB_ID"
printf 'source job     = %s\n' "$SOURCE_JOB_ID"
printf 'restart global = 100000\n'
printf 'segment length = 10000\n'
printf 'target global  = 110000\n'
printf 'wall model     = OFF\n'
printf '============================================================\n'

cd "$RUN_DIR"
set +e
srun --mpi=pmix --unbuffered --nodes=1 --ntasks=40 --ntasks-per-node=40 \
    --cpus-per-task=1 --distribution=block:block --cpu-bind=cores \
    --kill-on-bad-exit=1 stdbuf -oL -eL \
    "$EXE" flatplate "$STAGED" > solver.out 2> solver.err
RC=$?
set -e

printf '\n===== FINAL LOCAL ITERATION =====\n'
grep '^Iteration ' solver.out | tail -n 1 || true
printf '\n===== STDERR =====\n'
cat solver.err || true
printf '\nRETURN CODE = %s\n' "$RC"

[[ "$RC" -eq 0 ]] || { echo "CSAFM 0.50 SOLVER FAILED" >&2; exit "$RC"; }
if grep -Eqi 'nan|FATAL|ERROR:' solver.out solver.err; then
    echo "FATAL: non-finite/fatal diagnostic detected" >&2
    exit 3
fi
if grep -q 'SA production wall treatment: ON' solver.out; then
    echo "FATAL: wall treatment active" >&2
    exit 4
fi
grep -q '^Run complete$' solver.out || { echo "FATAL: Run complete missing" >&2; exit 5; }
grep -q '^Iteration 10000/10000 ' solver.out || { echo "FATAL: segment did not reach local 10000/10000" >&2; exit 6; }

FINAL_PVTU="$RUN_DIR/output/flatplate/flatplate_step_00110000.pvtu"
[[ -f "$FINAL_PVTU" ]] || { echo "FATAL: global-110000 PVTU missing: $FINAL_PVTU" >&2; exit 7; }
piece_count=$(find "$RUN_DIR/output/flatplate" -maxdepth 1 -type f -name 'flatplate_step_00110000_rank_*.vtu' | wc -l)
[[ "$piece_count" -eq 40 ]] || { echo "FATAL: expected 40 global-110000 VTU pieces; found $piece_count" >&2; exit 8; }
printf 'FINAL_PVTU=%s\n' "$FINAL_PVTU"

printf '\n===== TMR QUANTITATIVE DIAGNOSTICS =====\n'
python3 "$ANALYZER" "$STAGED" "$RUN_DIR/output/flatplate" \
    --station 0.97008 --re-x1 5.0e6 --u-inf 1.0 --kappa 0.41 \
    | tee "$JOB_ROOT/tmr_diagnostics.json"

printf '\nCSAFM_RESULT_ROOT=%s\n' "$JOB_ROOT"
printf 'CSAFM_DIAGNOSTICS=%s/tmr_diagnostics.json\n' "$JOB_ROOT"
printf 'CSAFM 0.50 SINGLE CASE: COMPLETED\n'
SLURM

chmod 700 "$JOB_SCRIPT"

printf '\n===== SUBMIT CSAFM 0.50 ONLY =====\n'
JOB_ID=$(sbatch --parsable \
    --export=ALL,SOURCE_JOB_ID="$SOURCE_JOB_ID",RESTART_ROOT="$RESTART_DIR" \
    "$JOB_SCRIPT")

printf 'JOBID         = %s\n' "$JOB_ID"
printf 'SOURCE JOB    = %s\n' "$SOURCE_JOB_ID"
printf 'CSAFM         = 0.50\n'
printf 'TARGET GLOBAL = 110000\n'
printf 'OUT           = %s/fp_c050_one_%s.out\n' "$LOGS" "$JOB_ID"
printf 'ERR           = %s/fp_c050_one_%s.err\n' "$LOGS" "$JOB_ID"
printf 'RUN ROOT      = %s/runs/sa_tmr_csafm_single/csafm_050_job_%s\n' "$CASE" "$JOB_ID"
printf '\nQueue:\n'
squeue -j "$JOB_ID" || true
printf '\nLive solver output after startup:\n'
printf 'tail -F %s/runs/sa_tmr_csafm_single/csafm_050_job_%s/run/solver.out\n' "$CASE" "$JOB_ID"
