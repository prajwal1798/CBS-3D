#!/bin/bash --login
#=============================================================================
# CBS3D++_SI OpenMP scaling study on Sunbird (Supercomputing Wales)
#
# Grabs ONE full compute node (40 cores, exclusive), then runs the same fixed
# workload at increasing OpenMP thread counts and records wall time for each.
# A fair strong-scaling test: identical work every run (fixed iteration count,
# early-stop disabled, VTU output disabled so we time compute, not disk).
#
# EDIT the two lines marked  <-- EDIT  then:   sbatch scaling_study.sh
# Results -> scaling_results.csv   (plot with plot_scaling.py)
#=============================================================================
#SBATCH --job-name=cbs3d_scaling
#SBATCH --output=scaling.out.%J
#SBATCH --error=scaling.err.%J
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --time=0-03:00
#SBATCH --account=scwXXXX                 # <-- EDIT: your SCW project code

module purge
module load compiler/gnu/9                # <-- EDIT: the GCC you build with (>=10 ideal for C++20)

set -euo pipefail

CASE=LidDrivenCavity3D
ITERS=300                                  # fixed iterations per run (fair timing)
THREAD_LIST="1 2 4 8 16 20 32 40"          # 20 = one socket, 40 = full node
ROOT=$(pwd)
BIN="$ROOT/cbs3dpp_si"

# OpenMP affinity (works on Sunbird/Linux; pins threads to cores for cache locality)
export OMP_PROC_BIND=close
export OMP_PLACES=cores

#-----------------------------------------------------------------------------
# 1. Build (uses your Makefile; -fopenmp already wired in)
#-----------------------------------------------------------------------------
echo "[build] make clean && make"
make clean >/dev/null
make

#-----------------------------------------------------------------------------
# 2. Build a TIMING input set: fixed iters, no early-stop, no VTU I/O
#    Edits (line:field) on a copy of the real .par:
#      line 10 f1  ntime              -> $ITERS
#      line 26 f1  vel_check          -> 0   (disable steady-state early stop)
#      line 44 f5  vtu_output_enabled -> 0   (no VTU writes during timing)
#      line 44 f8  write_boundary_dbg -> 0
#-----------------------------------------------------------------------------
WORK="$ROOT/scaling_run"
rm -rf "$WORK"; mkdir -p "$WORK"
cp "input/$CASE.plt" "input/$CASE.bco" "$WORK/"

awk -v iters="$ITERS" '
  NR==10 { $1=iters;             print; next }   # ntime
  NR==26 { $1=0;                 print; next }   # vel_check off
  NR==44 { $5=0; $8=0;           print; next }   # vtu off, debug off
  { print }
' "input/$CASE.par" > "$WORK/$CASE.par"

echo "[timing par] ntime line:"; sed -n '10p' "$WORK/$CASE.par"

#-----------------------------------------------------------------------------
# 3. Sweep thread counts, time each run
#-----------------------------------------------------------------------------
RESULTS="$ROOT/scaling_results.csv"
echo "threads,seconds,iters,sec_per_iter" > "$RESULTS"

cd "$WORK"
for NT in $THREAD_LIST; do
    export OMP_NUM_THREADS=$NT
    echo "[run] OMP_NUM_THREADS=$NT ..."
    t0=$(date +%s.%N)
    "$BIN" "$CASE" > "run_${NT}.log" 2>&1
    t1=$(date +%s.%N)
    dt=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')
    spi=$(awk -v d="$dt" -v n="$ITERS" 'BEGIN{printf "%.5f", d/n}')
    echo "$NT,$dt,$ITERS,$spi" >> "$RESULTS"
    echo "      -> ${dt}s  (${spi}s/iter)"
done

echo
echo "[done] results -> $RESULTS"
column -s, -t "$RESULTS"
