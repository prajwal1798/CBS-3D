#!/bin/bash -l
#SBATCH --job-name=cbs_p40_showcase
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=00:45:00
#SBATCH --output=logs/p40_showcase_%j.out
#SBATCH --error=logs/p40_showcase_%j.err

set -euo pipefail

REPO="/scratch/s.2337862/CBS3D++_SI_DD"
SUITE_ROOT="/scratch/s.2337862/CBS3D_cases/blanket/partition_sets/mesh_3eaed573_20260724"
PARTITION_ROOT="$SUITE_ROOT/P0040/cbs_partitions_40"
RUN_ROOT="/scratch/s.2337862/CBS3D_runs/blanket_p40_showcase_10000_20260726"
CASE_NAME="blanket"
RANKS=40
ITERATIONS=10000
EXE="$REPO/build-sunbird-dd/cbs3dpp_si"
BASE_TEMPLATE="$SUITE_ROOT/manifests/blanket_parameters_steady_10000.par"
JOB_ROOT="$RUN_ROOT/job_${SLURM_JOB_ID}"
STAGED_PARTITION_ROOT="$JOB_ROOT/partitions_40"
RUN_DIR="$JOB_ROOT/run"

mkdir -p "$JOB_ROOT" "$RUN_DIR"
source "$REPO/sunbird_modules.sh"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close
unset CBS3D_VERBOSE
unset CBS3D_ALL_RANK_OUTPUT

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ -x "$EXE" ]] || fail "missing executable: $EXE"
[[ -d "$PARTITION_ROOT" ]] || fail "missing P40 partition root: $PARTITION_ROOT"
[[ -f "$BASE_TEMPLATE" ]] || fail "missing parameter template: $BASE_TEMPLATE"
[[ "$SLURM_NTASKS" -eq "$RANKS" ]] || fail "expected 40 MPI ranks"
[[ "$SLURM_JOB_NUM_NODES" -eq 1 ]] || fail "expected one compute node"

SHOWCASE_TEMPLATE="$JOB_ROOT/blanket_steady_10000_vtu_every_1000.par"
python3 - "$BASE_TEMPLATE" "$SHOWCASE_TEMPLATE" <<'PY'
from pathlib import Path
import sys

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
lines = src.read_text(encoding="utf-8").splitlines()
if len(lines) < 44:
    raise SystemExit("parameter template is shorter than expected")
if "ntime, transient_on" not in lines[8]:
    raise SystemExit("unexpected parameter layout at line 9")
if "CBS3D output/monitor controls" not in lines[42]:
    raise SystemExit("unexpected parameter layout at line 43")

lines[9] = "10000  0  1  1.0e-5  1"
lines[43] = "1  100  100  0  1  1000  0.0  0  10000"
dst.write_bytes(("\n".join(lines) + "\n").encode("utf-8"))
PY

rm -rf "$STAGED_PARTITION_ROOT"
mkdir -p "$STAGED_PARTITION_ROOT"
cp -al "$PARTITION_ROOT/." "$STAGED_PARTITION_ROOT/"

PAR_COUNT=$(find "$STAGED_PARTITION_ROOT" -type f -name '*.par' | wc -l)
[[ "$PAR_COUNT" -eq "$RANKS" ]] || fail "found $PAR_COUNT parameter files, expected 40"

while IFS= read -r -d '' par_file; do
    cp --remove-destination "$SHOWCASE_TEMPLATE" "$par_file"
done < <(find "$STAGED_PARTITION_ROOT" -type f -name '*.par' -print0)

cd "$RUN_DIR"
start_epoch=$(date +%s.%N)

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
    "$EXE" "$CASE_NAME" "$STAGED_PARTITION_ROOT" \
    > solver.out \
    2> solver.err

end_epoch=$(date +%s.%N)
elapsed_time=$(awk -v a="$start_epoch" -v b="$end_epoch" 'BEGIN{printf "%.6f", b-a}')

OUTPUT_DIR="$RUN_DIR/output/$CASE_NAME"
[[ -d "$OUTPUT_DIR" ]] || fail "solver output directory was not created: $OUTPUT_DIR"

PVTU_COUNT=$(find "$OUTPUT_DIR" -maxdepth 1 -type f -name '*.pvtu' | wc -l)
VTU_COUNT=$(find "$OUTPUT_DIR" -maxdepth 1 -type f -name '*.vtu' | wc -l)
PVD_COUNT=$(find "$OUTPUT_DIR" -maxdepth 1 -type f -name '*.pvd' | wc -l)

cat > "$JOB_ROOT/summary.txt" <<SUMMARY
job_id=$SLURM_JOB_ID
partition=$SLURM_JOB_PARTITION
nodes=1
mpi_ranks=40
omp_threads=1
iterations=10000
vtu_output_enabled=1
vtu_output_every_iterations=1000
expected_snapshot_iterations=1000,2000,3000,4000,5000,6000,7000,8000,9000,10000
elapsed_time_s=$elapsed_time
pvtu_files=$PVTU_COUNT
vtu_files=$VTU_COUNT
pvd_files=$PVD_COUNT
output_directory=$OUTPUT_DIR
SUMMARY

[[ "$PVTU_COUNT" -ge 10 ]] || fail "expected at least 10 PVTU snapshots, found $PVTU_COUNT"
[[ "$VTU_COUNT" -ge 400 ]] || fail "expected at least 400 rank VTU files, found $VTU_COUNT"

echo "CBS3D P40 showcase run: PASS"
cat "$JOB_ROOT/summary.txt"
