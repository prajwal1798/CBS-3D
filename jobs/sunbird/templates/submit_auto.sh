#!/bin/bash
# ================================
# FILE      : submit_auto.sh
# PURPOSE   : LOOK AT WHAT IS FREE ON THE CLUSTER RIGHT NOW AND START THE
#             RUN ON THE LARGEST READY PARTITION SET THAT FITS:
#               1. LIST THE READY SETS: CASE_DIR/partitions_N (COMPLETE ONES)
#               2. COUNT USABLE IDLE PHYSICAL CORES PER NODE (IDLE CORES,
#                  LIMITED BY FREE MEMORY AT MEM_PER_RANK_GB PER RANK)
#               3. PICK THE LARGEST N THAT FITS (OR THE LARGEST <= --max)
#               4. WHOLE NODES (submit_mpi.sh ... nodes) WHEN N = k x 128 AND k
#                  NODES ARE COMPLETELY IDLE, OTHERWISE SCATTERED CORES (... any)
#
# USAGE     : bash submit_auto.sh CASE_DIR CASE [OPTIONS]
#   --time T      SLURM TIME LIMIT          (DEFAULT 24:00:00)
#   --max N       NEVER USE MORE THAN N RANKS
#   --dry-run     ONLY SHOW THE DECISION
#
# NOTES     : A SNAPSHOT: OTHER USERS MAY TAKE CORES BEFORE THE JOB STARTS;
#             THEN IT WAITS IN THE QUEUE. A RESTART MUST USE THE SAME N AS THE
#             CHECKPOINT, SO PASS --max/--dry-run AND submit_mpi.sh FOR THAT.
# ================================

set -euo pipefail
[[ "$#" -ge 2 ]] || { sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

source "$HOME/cbs3d_env.sh"
CASE_DIR="$(readlink -f "$1")"
CASE="$2"
shift 2
TIME="24:00:00"
MAX=0
DRY=0
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --time)    TIME="$2"; shift 2 ;;
        --max)     MAX="$2";  shift 2 ;;
        --dry-run) DRY=1;     shift ;;
        *) echo "ERROR: UNKNOWN OPTION $1" >&2; exit 1 ;;
    esac
done

T="$CBS3D_CODE/jobs/sunbird/templates"
PART_NAME="${CBS3D_PARTITION:-compute}"
CORES_PER_NODE="${CBS3D_CORES_PER_NODE:-128}"
TPC="${CBS3D_THREADS_PER_CORE:-2}"
MEM_PER_RANK_GB="${CBS3D_MEM_PER_RANK_GB:-4}"    # 2G PER CPU x 2 HYPERTHREADS

# ================================
# STEP 1: READY PARTITION SETS
# ================================
READY=()
for d in "$CASE_DIR"/partitions_*; do
    [[ -d "$d" ]] || continue
    N="${d##*_}"
    [[ "$N" =~ ^[0-9]+$ ]] || continue
    HAVE=$(find "$d" -maxdepth 1 -type d -name 'rank_[0-9][0-9][0-9][0-9]' | wc -l)
    [[ "$HAVE" -eq "$N" ]] && READY+=("$N")
done
[[ "${#READY[@]}" -gt 0 ]] || { echo "ERROR: NO COMPLETE partitions_N IN $CASE_DIR" >&2; exit 1; }
IFS=$'\n' READY=($(sort -n <<< "${READY[*]}")); unset IFS
echo "READY  : partitions ${READY[*]}"

# ================================
# STEP 2: USABLE IDLE CORES NOW
# ================================
read -r FREE WHOLE < <(
    sinfo -N -h -p "$PART_NAME" -O "NodeList:16,CPUsState:24,AllocMem:12,Memory:12" |
    awk -v tpc="$TPC" -v mpr="$MEM_PER_RANK_GB" '
    {
        split($2, c, "/")
        idle = int(c[2] / tpc)
        bymem = int((($4 - $3) / 1024) / mpr)
        use = (idle < bymem) ? idle : bymem
        free += use
        if (idle == int(c[4] / tpc) && use == idle && idle > 0) whole++
    }
    END { printf "%d %d\n", free, whole }')
echo "FREE   : $FREE usable idle cores, $WHOLE completely idle nodes ($PART_NAME)"

# ================================
# STEP 3: LARGEST SET THAT FITS
# ================================
PICK=0
for N in "${READY[@]}"; do
    [[ "$MAX" -gt 0 && "$N" -gt "$MAX" ]] && continue
    [[ "$N" -le "$FREE" ]] && PICK="$N"
done
if [[ "$PICK" -eq 0 ]]; then
    PICK="${READY[0]}"
    echo "NOTE   : NOT EVEN ${PICK} CORES FREE NOW; SUBMITTING ${PICK} RANKS TO WAIT IN THE QUEUE"
fi

NODES=$(( (PICK + CORES_PER_NODE - 1) / CORES_PER_NODE ))
# WHOLE NODES ONLY WHEN N FILLS THEM (N = k x 128) AND k NODES ARE IDLE;
# OTHERWISE SCATTERED CORES, SO NO EXCLUSIVE NODE SITS HALF EMPTY.
if [[ $(( PICK % CORES_PER_NODE )) -eq 0 && "$NODES" -le "$WHOLE" ]]; then MODE=nodes; else MODE=any; fi
echo "CHOICE : $PICK RANKS, MODE $MODE$( [[ $MODE == nodes ]] && echo " ($NODES WHOLE NODE(S))"), TIME $TIME"

[[ "$DRY" -eq 1 ]] && { echo "DRY RUN: NOTHING SUBMITTED"; exit 0; }
bash "$T/submit_mpi.sh" "$CASE_DIR" "$CASE" "$PICK" "$MODE" "$TIME"
