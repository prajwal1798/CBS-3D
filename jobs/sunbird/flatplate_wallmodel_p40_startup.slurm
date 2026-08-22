#!/bin/bash -l
#SBATCH --job-name=fp_wm_p40_start
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/s.2337862/CBS3D/cases/flatplate/jobs/logs/wallmodel_start_%j.out
#SBATCH --error=/scratch/s.2337862/CBS3D/cases/flatplate/jobs/logs/wallmodel_start_%j.err

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate

SOURCE_PART=$CASE/partition/coarse/P0040/cbs
EXE=$CODE/build/mpi-release/cbs3d_parallel

RUN_ROOT=$CASE/runs/coarse_wallmodel_startup
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

# Production momentum wall treatment.  Tangential velocity is no longer
# strongly clamped on BC530; impermeability remains strong and the Spalding
# shear replaces the natural viscous face flux in Step 1.
export CBS3D_SA_WALL_TREATMENT=1
export CBS3D_WALL_KAPPA=0.41
export CBS3D_WALL_B=5.2

# Startup diagnostics on rank 0 only.
export CBS3D_VERBOSE=1
unset CBS3D_ALL_RANK_OUTPUT

[[ -x "$EXE" ]] || {
    echo "ERROR: executable missing: $EXE" >&2
    exit 1
}

[[ -d "$SOURCE_PART" ]] || {
    echo "ERROR: partition root missing: $SOURCE_PART" >&2
    exit 1
}

[[ "$SLURM_NTASKS" -eq 40 ]] || {
    echo "ERROR: expected 40 MPI ranks" >&2
    exit 1
}

mkdir -p "$CASE/jobs/logs" "$JOB_ROOT" "$RUN_DIR"
rm -rf "$STAGED"
mkdir -p "$STAGED"

# Preserve the validated retained partition.  Every run gets a private staged
# copy so the source P0040 case remains immutable.
cp -a "$SOURCE_PART/." "$STAGED/"

rank_dirs=$(find "$STAGED" -maxdepth 1 -type d -name 'rank_[0-9][0-9][0-9][0-9]' | wc -l)
[[ "$rank_dirs" -eq 40 ]] || {
    echo "ERROR: expected 40 staged rank directories, found $rank_dirs" >&2
    exit 1
}

# Modify ONLY the staged parameter files for a ten-iteration P40 smoke run.
python3 - "$STAGED" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
files = sorted(root.glob("rank_*/*.par"))

if len(files) != 40:
    raise SystemExit(
        "ERROR: expected 40 .par files, found {}".format(len(files))
    )

for path in files:
    lines = path.read_text(encoding="utf-8").splitlines()

    got_time = False
    got_output = False
    got_sa = False

    for i, line in enumerate(lines):
        key = line.strip()

        if key == "ntime transient_on dtfixed dtfix iwrite":
            # Ten steady pseudo-time iterations.  Keep dtfixed=0 so the
            # retained flat-plate LTS controls remain active.
            lines[i + 1] = "10  0  0  1.0e-5  1"
            got_time = True

        if line.startswith("CBS3D output/monitor controls:"):
            # residual every iteration; console every iteration; VTU at 5/10;
            # boundary-debug arrays on.
            lines[i + 1] = "1  1  1  0  1  5  0.0  1  100"
            got_output = True

        if line.startswith("Spalart-Allmaras controls:"):
            fields = lines[i + 1].split()
            if len(fields) < 3 or fields[0] != "1" or fields[1] != "0":
                raise SystemExit(
                    "ERROR: expected standard SA enabled in {}".format(path)
                )
            got_sa = True

    if not got_time:
        raise SystemExit("ERROR: ntime block missing in {}".format(path))
    if not got_output:
        raise SystemExit("ERROR: monitor block missing in {}".format(path))
    if not got_sa:
        raise SystemExit("ERROR: SA control block missing in {}".format(path))

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

print("Prepared {} staged wall-model startup parameter files.".format(len(files)))
PY

printf '============================================================\n'
printf 'CBS3D FLAT-PLATE P40 WALL-MODEL STARTUP\n'
printf '============================================================\n'
printf 'date            = %s\n' "$(date -Is)"
printf 'host            = %s\n' "$(hostname)"
printf 'job             = %s\n' "$SLURM_JOB_ID"
printf 'MPI ranks       = %s\n' "$SLURM_NTASKS"
printf 'executable      = %s\n' "$EXE"
printf 'partition       = %s\n' "$STAGED"
printf 'run directory   = %s\n' "$RUN_DIR"
printf 'wall treatment  = %s\n' "$CBS3D_SA_WALL_TREATMENT"
printf 'kappa / B       = %s / %s\n' "$CBS3D_WALL_KAPPA" "$CBS3D_WALL_B"
printf 'git HEAD        = %s\n' "$(git -C "$CODE" rev-parse HEAD)"
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
echo "===== WALL-MODEL / TURBULENCE EVIDENCE ====="
grep -E \
    'SA production wall treatment|Spalding|RelSA|NuTilde|NuT=|MuT|MuEff|DivRMS|CG=|ERROR|FATAL|Run complete' \
    solver.out | tail -n 160 || true

echo
echo "===== LAST 160 SOLVER LINES ====="
tail -n 160 solver.out || true

echo
echo "===== STDERR ====="
cat solver.err || true

echo
echo "RETURN CODE = $RC"

if [[ "$RC" -ne 0 ]]; then
    echo "CBS3D FLAT-PLATE P40 WALL-MODEL STARTUP: FAILED" >&2
    exit "$RC"
fi

if ! grep -q 'SA production wall treatment: ON' solver.out; then
    echo "ERROR: solver completed without confirming production wall treatment." >&2
    exit 2
fi

if grep -Eqi 'NaN|nan|Inf|FATAL|ERROR:' solver.out solver.err; then
    echo "ERROR: non-finite/fatal diagnostic detected in startup logs." >&2
    exit 3
fi

echo "CBS3D FLAT-PLATE P40 WALL-MODEL STARTUP: COMPLETED"
