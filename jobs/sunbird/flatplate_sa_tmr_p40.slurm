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

# Mandatory wall-resolved benchmark: no wall model may leak into this run.
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

# Patch the retained .par files using the SAME positional grammar as MeshIO.
# The old version searched for human-readable labels.  The retained partition
# predates those exact label strings even though the C++ parser reads it just
# fine positionally, which is why job 8197461 died before srun.
python3 - "$STAGED" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
files = sorted(root.glob("rank_*/*.par"))
if len(files) != 40:
    raise SystemExit("FATAL: expected 40 .par files, found {}".format(len(files)))

def active_indices(lines):
    return [
        i for i, line in enumerate(lines)
        if line.strip() and not line.lstrip().startswith(("#", "!"))
    ]

def find_active_position(lines, active, predicate, description):
    hits = [p for p, i in enumerate(active) if predicate(lines[i].strip().lower())]
    if len(hits) != 1:
        raise SystemExit(
            "FATAL: expected one {} anchor; found {}".format(description, len(hits))
        )
    return hits[0]

for path in files:
    lines = path.read_text(encoding="utf-8").splitlines()
    active = active_indices(lines)

    # MeshIO reads label/data pairs.  Anchor on labels that are known to exist
    # in the retained files, then locate legacy-labelled fields by their parser
    # position instead of their prose label.
    time_pos = find_active_position(
        lines, active,
        lambda s: s == "ntime transient_on dtfixed dtfix iwrite",
        "ntime")
    out_pos = find_active_position(
        lines, active,
        lambda s: s.startswith("cbs3d output/monitor controls:"),
        "output/monitor")

    # From MeshIO grammar:
    # ntime block -> csafm -> timestep -> Re -> convection -> pnode ->
    # write-output -> rel/abs tolerance -> steady controls.
    steady_data_pos = time_pos + 17

    # output/monitor block -> artificial-diffusion block -> SA control block.
    art_data_pos = out_pos + 3
    sa_data_pos = out_pos + 5

    needed = max(time_pos + 1, steady_data_pos, out_pos + 1, art_data_pos, sa_data_pos)
    if needed >= len(active):
        raise SystemExit("FATAL: positional .par grammar truncated in {}".format(path))

    time_i = active[time_pos + 1]
    steady_i = active[steady_data_pos]
    output_i = active[out_pos + 1]
    art_i = active[art_data_pos]
    sa_i = active[sa_data_pos]

    time_values = lines[time_i].split()
    if len(time_values) < 5:
        raise SystemExit("FATAL: malformed timestep row in {}".format(path))
    time_values[0] = "100000"
    time_values[1] = "0"
    time_values[4] = "1000"
    lines[time_i] = "  ".join(time_values)

    steady_values = lines[steady_i].split()
    if len(steady_values) < 4:
        raise SystemExit("FATAL: malformed steady-control row in {}".format(path))
    # Preserve the velocity/temperature controls and force SA convergence into
    # the stop rule: sa_check=1, RelSA tolerance=1e-7.
    lines[steady_i] = "  ".join(steady_values[:4] + ["1", "1.0e-7"])

    # residual every 10; console every 1000; VTU every 10000; no live plot;
    # boundary diagnostics retained; minimum 5000 iterations.
    lines[output_i] = "1  10  1000  0  1  10000  0.0  1  5000"

    art_values = lines[art_i].split()
    if not art_values or int(float(art_values[0])) != 0:
        raise SystemExit(
            "FATAL: wall-resolved benchmark requires art_diff=0 in {}".format(path)
        )

    sa_values = lines[sa_i].split()
    if len(sa_values) < 3 or sa_values[:3] != ["1", "0", "0"]:
        raise SystemExit(
            "FATAL: expected standard SA [1 0 0] in {}; got {}".format(path, sa_values)
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

print("Prepared 40 staged wall-resolved TMR parameter files: PASS")
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
if [[ -r "$GIT_HEAD_FILE" ]]; then GIT_HEAD=$(cat "$GIT_HEAD_FILE"); else GIT_HEAD="unavailable-on-compute-node"; fi
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
grep -E 'Iter|RelU|RelV|RelW|RelP|RelSA|NuTilde|NuT=|MuT|MuEff|DivRMS|CG=|Run complete|FATAL|ERROR' solver.out | tail -n 240 || true

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
[[ -n "$FINAL_PVTU" ]] || { echo "FATAL: no distributed VTU snapshot found." >&2; exit 5; }

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
