#!/bin/bash -l
#SBATCH --job-name=fp_sa_tmr_p40
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/s.2337862/CBS3D/cases/flatplate/jobs/logs/sa_tmr_%j.out
#SBATCH --error=/scratch/s.2337862/CBS3D/cases/flatplate/jobs/logs/sa_tmr_%j.err

# ==============================================================================
# Wall-resolved Spalart-Allmaras NASA/TMR flat-plate benchmark.
#
# This job deliberately does NOT enable either momentum wall model.  It is the
# numerical-verification gate for the SA transport + CBS/RANS coupling itself.
# The retained P0040 partition is copied to a job-private staging tree before
# parameter changes, so the source partition remains immutable.
# ==============================================================================

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
SOURCE_PART=$CASE/partition/coarse/P0040/cbs
EXE=$CODE/build/mpi-release/cbs3d_parallel
ANALYZER=$CODE/examples/Flat_Plate_Turbulent/analyze_distributed_tmr.py
RUN_ROOT=$CASE/runs/sa_tmr_wall_resolved
JOB_ROOT=$RUN_ROOT/job_${SLURM_JOB_ID}
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

# This is the wall-resolved benchmark.  If either wall law leaks into this job,
# the comparison no longer tests the SA finite-element discretisation.
unset CBS3D_SA_WALL_TREATMENT
unset CBS3D_CHT_WALL_TREATMENT
unset CBS3D_WALL_KAPPA
unset CBS3D_WALL_B

export CBS3D_VERBOSE=1
unset CBS3D_ALL_RANK_OUTPUT

[[ -x "$EXE" ]] || { echo "FATAL: executable missing: $EXE" >&2; exit 1; }
[[ -f "$ANALYZER" ]] || { echo "FATAL: TMR analyzer missing: $ANALYZER" >&2; exit 1; }
[[ -d "$SOURCE_PART" ]] || { echo "FATAL: partition missing: $SOURCE_PART" >&2; exit 1; }
[[ "$SLURM_NTASKS" -eq 40 ]] || { echo "FATAL: expected 40 MPI ranks" >&2; exit 1; }

mkdir -p "$CASE/jobs/logs" "$JOB_ROOT" "$RUN_DIR"
rm -rf "$STAGED"
mkdir -p "$STAGED"
cp -a "$SOURCE_PART/." "$STAGED/"

rank_dirs=$(find "$STAGED" -maxdepth 1 -type d -name 'rank_[0-9][0-9][0-9][0-9]' | wc -l)
[[ "$rank_dirs" -eq 40 ]] || { echo "FATAL: expected 40 staged rank directories, found $rank_dirs" >&2; exit 1; }

# Max 100k iterations for direct comparison with the previous benchmark.  The
# run may terminate earlier, but only after both velocity and SA residual gates
# have passed.  Output cadence is intentionally sparse enough not to turn this
# CFD benchmark into an I/O benchmark, one of computing's less glorious genres.
python3 - "$STAGED" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
files = sorted(root.glob("rank_*/*.par"))
if len(files) != 40:
    raise SystemExit("FATAL: expected 40 .par files, found {}".format(len(files)))

for path in files:
    lines = path.read_text(encoding="utf-8").splitlines()
    found_time = found_steady = found_output = found_sa = found_art = False

    for i, line in enumerate(lines):
        key = line.strip()
        lower = key.lower()

        if key == "ntime transient_on dtfixed dtfix iwrite":
            values = lines[i + 1].split()
            if len(values) < 5:
                raise SystemExit("FATAL: malformed timestep row in {}".format(path))
            values[0] = "100000"
            values[1] = "0"
            # Preserve retained LTS/fixed-dt selection and dt value; only cap the
            # benchmark length and reduce legacy output chatter.
            values[4] = "1000"
            lines[i + 1] = "  ".join(values)
            found_time = True

        elif "velocity/temperature steady controls" in lower:
            values = lines[i + 1].split()
            if len(values) < 4:
                raise SystemExit("FATAL: malformed steady-control row in {}".format(path))
            # Preserve velocity/temperature controls, but make SA convergence a
            # mandatory steady-stop condition.
            lines[i + 1] = "  ".join(values[:4] + ["1", "1.0e-7"])
            found_steady = True

        elif key.startswith("CBS3D output/monitor controls:"):
            # residual CSV every 10 iters; console every 1000; VTU every 10000;
            # no live plot; retain boundary-debug fields; no early stop before
            # 5000 iterations even if the initial field happens to look quiet.
            lines[i + 1] = "1  10  1000  0  1  10000  0.0  1  5000"
            found_output = True

        elif key.startswith("artificial diffusion control:"):
            values = lines[i + 1].split()
            if not values or values[0] != "0":
                raise SystemExit("FATAL: wall-resolved benchmark requires art_diff=0 in {}".format(path))
            found_art = True

        elif key.startswith("Spalart-Allmaras controls:"):
            values = lines[i + 1].split()
            if len(values) < 3 or values[:3] != ["1", "0", "0"]:
                raise SystemExit(
                    "FATAL: expected standard SA [1 0 0] in {}; got {}".format(path, values)
                )
            found_sa = True

    missing = [
        name for name, ok in (
            ("ntime", found_time),
            ("steady controls", found_steady),
            ("output controls", found_output),
            ("artificial diffusion", found_art),
            ("SA controls", found_sa),
        ) if not ok
    ]
    if missing:
        raise SystemExit("FATAL: {} missing {}".format(path, ", ".join(missing)))

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

print("Prepared 40 staged wall-resolved TMR parameter files.")
PY

printf '============================================================\n'
printf 'CBS3D WALL-RESOLVED SA NASA/TMR FLAT-PLATE BENCHMARK\n'
printf '============================================================\n'
printf 'date            = %s\n' "$(date -Is)"
printf 'host            = %s\n' "$(hostname)"
printf 'job             = %s\n' "$SLURM_JOB_ID"
printf 'MPI ranks       = %s\n' "$SLURM_NTASKS"
printf 'executable      = %s\n' "$EXE"
printf 'partition       = %s\n' "$STAGED"
printf 'run directory   = %s\n' "$RUN_DIR"
printf 'wall model      = OFF (mandatory)\n'
printf 'max iterations  = 100000\n'
printf 'SA stop gate    = 1.0e-7 after >=5000 iterations\n'
GIT_HEAD_FILE="$CODE/.git/refs/heads/integration/sa-mpi"
if [[ -r "$GIT_HEAD_FILE" ]]; then
    GIT_HEAD=$(cat "$GIT_HEAD_FILE")
else
    GIT_HEAD="unavailable-on-compute-node"
fi
printf 'git HEAD        = %s\n' "$GIT_HEAD"
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
    'Iter|RelU|RelV|RelW|RelP|RelSA|NuTilde|NuT=|MuT|MuEff|DivRMS|CG=|Run complete|FATAL|ERROR' \
    solver.out | tail -n 240 || true

echo
echo "===== LAST 160 SOLVER LINES ====="
tail -n 160 solver.out || true

echo
echo "===== STDERR ====="
cat solver.err || true

echo
echo "RETURN CODE = $RC"

if [[ "$RC" -ne 0 ]]; then
    echo "CBS3D WALL-RESOLVED SA TMR BENCHMARK: FAILED" >&2
    exit "$RC"
fi

if grep -Eqi 'nan|FATAL|ERROR:' solver.out solver.err; then
    echo "FATAL: non-finite/fatal diagnostic detected in benchmark logs." >&2
    exit 3
fi

if grep -q 'SA production wall treatment: ON' solver.out; then
    echo "FATAL: Spalding wall treatment was active in a wall-resolved benchmark." >&2
    exit 4
fi

FINAL_PVTU=$(find "$RUN_DIR/output/flatplate" -maxdepth 1 -type f -name '*.pvtu' 2>/dev/null | sort | tail -1 || true)
printf 'FINAL_PVTU=%s\n' "$FINAL_PVTU"
printf 'BENCHMARK_ROOT=%s\n' "$JOB_ROOT"

if [[ -z "$FINAL_PVTU" ]]; then
    echo "FATAL: solver completed but no distributed VTU snapshot was found." >&2
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
echo "CBS3D WALL-RESOLVED SA TMR BENCHMARK: COMPLETED"
