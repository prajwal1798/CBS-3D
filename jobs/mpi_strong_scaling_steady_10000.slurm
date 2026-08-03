#!/bin/bash -l

set -euo pipefail

: "${CBS3D_RANKS:?CBS3D_RANKS is required}"
: "${CBS3D_NODES:?CBS3D_NODES is required}"
: "${CBS3D_PARTITION_ROOT:?CBS3D_PARTITION_ROOT is required}"
: "${CBS3D_RUN_ROOT:?CBS3D_RUN_ROOT is required}"

REPO="${CBS3D_REPO:-/scratch/s.2337862/CBS3D++_SI_DD}"
SUITE_ROOT="${CBS3D_SUITE_ROOT:-/scratch/s.2337862/CBS3D_cases/blanket/partition_sets/mesh_3eaed573_20260724}"
CASE_NAME="${CBS3D_CASE_NAME:-blanket}"
ITERATIONS="${CBS3D_ITERATIONS:-10000}"
REPEATS="${CBS3D_REPEATS:-4}"
RANKS="${CBS3D_RANKS}"
NODES="${CBS3D_NODES}"
RANKS_PER_NODE=40
EXE="$REPO/build-sunbird-dd/cbs3dpp_si"
BASE_TEMPLATE="$SUITE_ROOT/manifests/blanket_parameters_steady_10000.par"
JOB_ROOT="$CBS3D_RUN_ROOT/job_${SLURM_JOB_ID}"
STAGED_PARTITION_ROOT="$JOB_ROOT/partitions_${RANKS}"
TIMINGS_CSV="$JOB_ROOT/timings.csv"
PLACEMENT_LOG="$JOB_ROOT/placement.log"

mkdir -p "$JOB_ROOT"

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
[[ -d "$CBS3D_PARTITION_ROOT" ]] || fail "missing partition root: $CBS3D_PARTITION_ROOT"
[[ -f "$BASE_TEMPLATE" ]] || fail "missing parameter template: $BASE_TEMPLATE"
[[ "$SLURM_NTASKS" -eq "$RANKS" ]] || fail "SLURM_NTASKS=$SLURM_NTASKS, expected $RANKS"
[[ "$SLURM_JOB_NUM_NODES" -eq "$NODES" ]] || fail "allocated nodes=$SLURM_JOB_NUM_NODES, expected $NODES"
[[ $((NODES * RANKS_PER_NODE)) -eq "$RANKS" ]] || fail "expected 40 ranks per node"
[[ "$ITERATIONS" -eq 10000 ]] || fail "this campaign is fixed at 10000 iterations"
[[ "$REPEATS" -ge 2 ]] || fail "REPEATS must include one warm-up and at least one measured run"

BENCH_TEMPLATE="$JOB_ROOT/blanket_steady_10000_noio.par"
python3 - "$BASE_TEMPLATE" "$BENCH_TEMPLATE" <<'PY'
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

# Exact 10,000-iteration steady-state benchmark.
lines[9] = "10000  0  1  1.0e-5  1"
# Sparse residual and console logging; VTU and debug output disabled.
lines[43] = "1  1000  1000  0  0  10000  0.0  0  10000"
dst.write_bytes(("\n".join(lines) + "\n").encode("utf-8"))
PY

rm -rf "$STAGED_PARTITION_ROOT"
mkdir -p "$STAGED_PARTITION_ROOT"
cp -al "$CBS3D_PARTITION_ROOT/." "$STAGED_PARTITION_ROOT/"

PAR_COUNT=$(find "$STAGED_PARTITION_ROOT" -type f -name '*.par' | wc -l)
[[ "$PAR_COUNT" -eq "$RANKS" ]] || fail "found $PAR_COUNT rank parameter files, expected $RANKS"

while IFS= read -r -d '' par_file; do
    cp --remove-destination "$BENCH_TEMPLATE" "$par_file"
done < <(find "$STAGED_PARTITION_ROOT" -type f -name '*.par' -print0)

while IFS= read -r -d '' par_file; do
    cmp -s "$BENCH_TEMPLATE" "$par_file" || fail "parameter mismatch: $par_file"
done < <(find "$STAGED_PARTITION_ROOT" -type f -name '*.par' -print0)

srun \
    --mpi=pmix \
    --nodes="$NODES" \
    --ntasks="$RANKS" \
    --ntasks-per-node="$RANKS_PER_NODE" \
    --cpus-per-task=1 \
    --distribution=block:block \
    --cpu-bind=cores \
    --kill-on-bad-exit=1 \
    bash -lc 'printf "rank=%s node=%s local_rank=%s\n" "$SLURM_PROCID" "$(hostname)" "$SLURM_LOCALID"' \
    > "$PLACEMENT_LOG"

PLACEMENT_TASKS=$(wc -l < "$PLACEMENT_LOG")
PLACEMENT_NODES=$(awk '{for (i=1;i<=NF;i++) if ($i ~ /^node=/) {sub(/^node=/,"",$i); print $i}}' "$PLACEMENT_LOG" | sort -u | wc -l)
MIN_PER_NODE=$(awk '{for (i=1;i<=NF;i++) if ($i ~ /^node=/) {sub(/^node=/,"",$i); print $i}}' "$PLACEMENT_LOG" | sort | uniq -c | awk 'NR==1{m=$1} $1<m{m=$1} END{print m+0}')
MAX_PER_NODE=$(awk '{for (i=1;i<=NF;i++) if ($i ~ /^node=/) {sub(/^node=/,"",$i); print $i}}' "$PLACEMENT_LOG" | sort | uniq -c | awk 'NR==1{m=$1} $1>m{m=$1} END{print m+0}')

[[ "$PLACEMENT_TASKS" -eq "$RANKS" ]] || fail "placement task count=$PLACEMENT_TASKS, expected $RANKS"
[[ "$PLACEMENT_NODES" -eq "$NODES" ]] || fail "placement node count=$PLACEMENT_NODES, expected $NODES"
[[ "$MIN_PER_NODE" -eq 40 && "$MAX_PER_NODE" -eq 40 ]] || fail "placement is not 40 ranks per node"

cat > "$JOB_ROOT/metadata.txt" <<META
case=$CASE_NAME
job_id=$SLURM_JOB_ID
partition=$SLURM_JOB_PARTITION
nodes=$NODES
mpi_ranks=$RANKS
ranks_per_node=$RANKS_PER_NODE
omp_threads=1
iterations=$ITERATIONS
simulation_mode=steady_state
transient_on=0
warmup_runs=1
measured_runs=$((REPEATS - 1))
mesh_tetrahedra=2383581
mesh_nodes=407458
mesh_boundary_triangles=28328
mesh_sha256=3eaed57362bc68144d9fc0ff9bc1dd898e112e3df116d2bb18fad125e39e0136
executable_sha256=$(sha256sum "$EXE" | awk '{print $1}')
parameter_sha256=$(sha256sum "$BENCH_TEMPLATE" | awk '{print $1}')
partition_root=$CBS3D_PARTITION_ROOT
staged_partition_root=$STAGED_PARTITION_ROOT
slurm_nodelist=$SLURM_JOB_NODELIST
created=$(date --iso-8601=seconds)
META

printf 'run,classification,setup_time_s,loop_time_s,elapsed_time_s,sec_per_iteration\n' > "$TIMINGS_CSV"

for ((run_index=0; run_index<REPEATS; ++run_index)); do
    if [[ "$run_index" -eq 0 ]]; then
        classification="warmup"
    else
        classification="measured"
    fi

    RUN_DIR="$JOB_ROOT/run_$(printf '%02d' "$run_index")"
    mkdir -p "$RUN_DIR"
    rm -rf "$RUN_DIR/output"
    cd "$RUN_DIR"

    start_epoch=$(date +%s.%N)
    set +e
    srun \
        --mpi=pmix \
        --unbuffered \
        --nodes="$NODES" \
        --ntasks="$RANKS" \
        --ntasks-per-node="$RANKS_PER_NODE" \
        --cpus-per-task=1 \
        --distribution=block:block \
        --cpu-bind=cores \
        --kill-on-bad-exit=1 \
        stdbuf -oL -eL \
        "$EXE" "$CASE_NAME" "$STAGED_PARTITION_ROOT" \
        > solver.out \
        2> solver.err
    solver_rc=$?
    set -e
    end_epoch=$(date +%s.%N)

    [[ "$solver_rc" -eq 0 ]] || fail "run $run_index failed with exit code $solver_rc"

    setup_time=$(awk '/setup time/ {value=$(NF-1)} END{print value}' solver.out)
    loop_time=$(awk '/loop time/ {value=$(NF-1)} END{print value}' solver.out)
    elapsed_time=$(awk -v a="$start_epoch" -v b="$end_epoch" 'BEGIN{printf "%.6f", b-a}')

    [[ -n "$setup_time" ]] || fail "setup-time extraction failed for run $run_index"
    [[ -n "$loop_time" ]] || fail "loop-time extraction failed for run $run_index"

    sec_per_iteration=$(awk -v t="$loop_time" -v n="$ITERATIONS" 'BEGIN{printf "%.12f", t/n}')
    printf '%d,%s,%s,%s,%s,%s\n' \
        "$run_index" "$classification" "$setup_time" "$loop_time" "$elapsed_time" "$sec_per_iteration" \
        >> "$TIMINGS_CSV"
done

python3 - "$TIMINGS_CSV" "$JOB_ROOT/statistics.txt" <<'PY'
import csv
import statistics
import sys

csv_path, out_path = sys.argv[1:3]
with open(csv_path, newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))
measured = [float(row["loop_time_s"]) for row in rows if row["classification"] == "measured"]
if not measured:
    raise SystemExit("no measured runs found")
with open(out_path, "w", encoding="utf-8") as handle:
    handle.write(f"measured_runs={len(measured)}\n")
    handle.write(f"minimum_loop_time_s={min(measured):.9f}\n")
    handle.write(f"median_loop_time_s={statistics.median(measured):.9f}\n")
    handle.write(f"mean_loop_time_s={statistics.mean(measured):.9f}\n")
    handle.write(f"maximum_loop_time_s={max(measured):.9f}\n")
    handle.write(f"population_stddev_s={statistics.pstdev(measured):.9f}\n")
PY

echo "CBS3D MPI strong-scaling run: PASS"
echo "Results: $JOB_ROOT"
cat "$JOB_ROOT/statistics.txt"
