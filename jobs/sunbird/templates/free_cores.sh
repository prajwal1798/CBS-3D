#!/bin/bash
# ================================
# FILE      : free_cores.sh
# PURPOSE   : HOW MANY PHYSICAL CORES ARE IDLE RIGHT NOW, NODE BY NODE.
#             SLURM COUNTS HYPERTHREADS (256 PER SUNBIRD2 NODE); THIS
#             PRINTS PHYSICAL CORES (128 PER NODE) AND FREE MEMORY.
#
# USAGE     : bash free_cores.sh [PARTITION]      (DEFAULT: compute)
#
# NOTES     : A SNAPSHOT: OTHER JOBS MAY TAKE THE CORES BEFORE YOURS STARTS.
#             IDLE CORES ON A NODE WHOSE MEMORY IS FULL CANNOT BE USED.
# ================================

set -euo pipefail
PART="${1:-compute}"
TPC="${CBS3D_THREADS_PER_CORE:-2}"

sinfo -N -h -p "$PART" -O "NodeList:16,CPUsState:24,AllocMem:12,Memory:12,StateCompact:10" |
awk -v tpc="$TPC" '
{
    split($2, c, "/")                       # ALLOCATED/IDLE/OTHER/TOTAL CPUS
    idle = int(c[2] / tpc)
    free_gb = int(($4 - $3) / 1024)
    if (idle > 0)
    {
        printf "%-12s %4d idle cores  %5d GB free  %s\n", $1, idle, free_gb, $5
        total += idle
        if (idle == int(c[4] / tpc)) whole++
    }
}
END {
    printf "TOTAL: %d idle physical cores, %d completely idle nodes\n", total, whole
}'
