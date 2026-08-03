#!/bin/bash -l
#SBATCH --job-name=cbs3d-dd5a
#SBATCH --partition=compute
#SBATCH --time=00:15:00
#SBATCH --nodes=12
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=logs/dd5a_%j.out
#SBATCH --error=logs/dd5a_%j.err

set -euo pipefail

: "${DD5A_PARTITION_ROOT:?Submit with --export=ALL,DD5A_PARTITION_ROOT=/path/to/prepared/tree}"

CASE_NAME="${CASE_NAME:-blanket}"
EXPECTED_RANKS="${EXPECTED_RANKS:-40}"
EXPECTED_ITERATIONS="${EXPECTED_ITERATIONS:-100}"
REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"

cd "${REPOSITORY_ROOT}"
mkdir -p logs output

source ./sunbird_modules.sh

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

RUN_TAG="job_${SLURM_JOB_ID}"
PLACEMENT_LOG="logs/dd5a_${SLURM_JOB_ID}_placement.log"
RUN_LOG="logs/dd5a_${SLURM_JOB_ID}.out"
RESIDUAL_CSV="output/${CASE_NAME}/${CASE_NAME}_distributed_residuals.csv"

if [[ ! -x ./build-sunbird-dd/cbs3dpp_si ]]; then
    echo "ERROR: solver executable is missing: ./build-sunbird-dd/cbs3dpp_si" >&2
    exit 2
fi

if [[ ! -d "${DD5A_PARTITION_ROOT}" ]]; then
    echo "ERROR: prepared partition root is missing: ${DD5A_PARTITION_ROOT}" >&2
    exit 2
fi

if [[ -d "output/${CASE_NAME}" ]]; then
    ARCHIVE="output/${CASE_NAME}_before_${RUN_TAG}_$(date +%Y%m%d_%H%M%S)"
    mv "output/${CASE_NAME}" "${ARCHIVE}"
    echo "Archived pre-existing output: ${ARCHIVE}"
fi

echo "============================================================"
echo "CBS3D DD-5A BATCH VALIDATION"
echo "============================================================"
echo "Slurm job                  : ${SLURM_JOB_ID}"
echo "allocated nodes            : ${SLURM_JOB_NUM_NODES}"
echo "allocated task slots       : ${SLURM_NTASKS}"
echo "solver MPI ranks           : ${EXPECTED_RANKS}"
echo "requested iterations       : ${EXPECTED_ITERATIONS}"
echo "case                       : ${CASE_NAME}"
echo "partition root             : ${DD5A_PARTITION_ROOT}"
echo "============================================================"

# The batch allocation reserves four CPU slots on every node. The solver step
# launches forty ranks cyclically over those forty-eight slots, producing the
# required balanced 3--4 ranks per node while preserving the forty-partition
# communicator.
srun \
    --nodes=12 \
    --ntasks="${EXPECTED_RANKS}" \
    --cpus-per-task=1 \
    --cpu-bind=cores \
    --distribution=cyclic \
    --kill-on-bad-exit=1 \
    bash -lc \
    'printf "PLACEMENT rank=%s node=%s local_rank=%s\n" \
        "$SLURM_PROCID" "$(hostname)" "$SLURM_LOCALID"' \
    > "${PLACEMENT_LOG}"

ACTIVE_NODES="$(
    sed -E 's/.*node=([^ ]+).*/\1/' "${PLACEMENT_LOG}" \
        | sort -u \
        | wc -l
)"

ACTIVE_TASKS="$(grep -c '^PLACEMENT' "${PLACEMENT_LOG}")"

MIN_TASKS_PER_NODE="$(
    sed -E 's/.*node=([^ ]+).*/\1/' "${PLACEMENT_LOG}" \
        | sort \
        | uniq -c \
        | awk 'NR==1 {min=$1} $1<min {min=$1} END {print min}'
)"

MAX_TASKS_PER_NODE="$(
    sed -E 's/.*node=([^ ]+).*/\1/' "${PLACEMENT_LOG}" \
        | sort \
        | uniq -c \
        | awk 'NR==1 {max=$1} $1>max {max=$1} END {print max}'
)"

echo "Placement active nodes     : ${ACTIVE_NODES}"
echo "Placement active tasks     : ${ACTIVE_TASKS}"
echo "Placement minimum/node     : ${MIN_TASKS_PER_NODE}"
echo "Placement maximum/node     : ${MAX_TASKS_PER_NODE}"

if [[ "${ACTIVE_NODES}" -ne 12 ||
      "${ACTIVE_TASKS}" -ne "${EXPECTED_RANKS}" ||
      "${MIN_TASKS_PER_NODE}" -ne 3 ||
      "${MAX_TASKS_PER_NODE}" -ne 4 ]]; then
    echo "ERROR: DD-5A rank placement is not balanced across twelve nodes" >&2
    exit 3
fi

echo "DD-5A placement validation : PASS"

set +e
srun \
    --unbuffered \
    --nodes=12 \
    --ntasks="${EXPECTED_RANKS}" \
    --cpus-per-task=1 \
    --cpu-bind=cores \
    --distribution=cyclic \
    --mpi=pmix \
    --kill-on-bad-exit=1 \
    stdbuf -oL -eL \
    ./build-sunbird-dd/cbs3dpp_si \
    "${CASE_NAME}" \
    "${DD5A_PARTITION_ROOT}"
SOLVER_RC=$?
set -e

echo "Solver return code         : ${SOLVER_RC}"

if [[ "${SOLVER_RC}" -ne 0 ]]; then
    echo "ERROR: DD-5A solver step failed" >&2
    exit "${SOLVER_RC}"
fi

if [[ ! -f "${RESIDUAL_CSV}" ]]; then
    echo "ERROR: distributed residual CSV was not created: ${RESIDUAL_CSV}" >&2
    exit 4
fi

EXPECTED_LINES=$((EXPECTED_ITERATIONS + 1))
ACTUAL_LINES="$(wc -l < "${RESIDUAL_CSV}")"
FIELD_FILE_COUNT="$(
    find "output/${CASE_NAME}" \
        -maxdepth 1 \
        -type f \
        \( -name '*.vtu' -o -name '*.pvtu' -o -name '*.pvd' \) \
        | wc -l
)"

echo "Residual CSV lines         : ${ACTUAL_LINES}"
echo "Expected CSV lines         : ${EXPECTED_LINES}"
echo "VTU/PVTU/PVD file count    : ${FIELD_FILE_COUNT}"

if [[ "${ACTUAL_LINES}" -ne "${EXPECTED_LINES}" ]]; then
    echo "ERROR: DD-5A residual history is incomplete" >&2
    exit 5
fi

if [[ "${FIELD_FILE_COUNT}" -ne 0 ]]; then
    echo "ERROR: DD-5A unexpectedly produced field-output files" >&2
    exit 6
fi

# Flush Slurm's stdout file before the checker inspects the production-loop
# completion banner and persistent PETSc build counts.
sync
sleep 1

python3 \
    tools/dd5a_check_stability.py \
    "${RESIDUAL_CSV}" \
    "${RUN_LOG}" \
    --expected-iterations "${EXPECTED_ITERATIONS}" \
    --expected-ranks "${EXPECTED_RANKS}"

echo "============================================================"
echo "CBS3D DD-5A BATCH VALIDATION: PASS"
echo "Residual CSV              : ${RESIDUAL_CSV}"
echo "Placement log             : ${PLACEMENT_LOG}"
echo "Slurm stdout              : ${RUN_LOG}"
echo "Slurm stderr              : logs/dd5a_${SLURM_JOB_ID}.err"
echo "============================================================"
