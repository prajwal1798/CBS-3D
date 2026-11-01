#!/bin/bash -l
# ================================
# FILE      : selfcheck_slab.slurm
# PURPOSE   : END-TO-END CHECK OF THE BUILT SOLVERS ON THE nuclear_heating_slab
#             CASE (EXACT SOLUTION KNOWN):
#               1. MESH slab.geo, CONVERT TO A SERIAL CBS CASE
#               2. SERIAL RUN    (CBS3D_SERIAL_EXE)
#               3. 4-RANK MPI RUN (CBS3D_MPI_EXE ON 4 PARTITIONS)
#                  (WITH CBS3D_TEMPERATURE_AFC=1, AS IN PRODUCTION)
#               4. BOTH AGAINST THE EXACT 1D TEMPERATURE (tools/check_nuclear_slab.py),
#                  AND THE MPI FIELD AGAINST THE SERIAL FIELD, NODE BY NODE
#               5. THE PRINTED HEAT INPUT AGAINST ITS EXACT VALUE (16.92 W)
#             LAST LINE: SELF-CHECK PASS OR SELF-CHECK FAIL.
#
# USAGE     : sbatch selfcheck_slab.slurm          (selfcheck.sh SUBMITS IT)
#
# OUTPUT    : $CBS3D_SCRATCH/selfcheck/<JOBID>/ AND cbs3d_selfcheck_<JOBID>.out
# ================================
#SBATCH --job-name=cbs3d_selfcheck
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=cbs3d_selfcheck_%j.out
#SBATCH --error=cbs3d_selfcheck_%j.err

set -uo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1

SRC="$CBS3D_CODE/cases/nuclear_heating_slab"
WORK="$CBS3D_SCRATCH/selfcheck/${SLURM_JOB_ID:-local}"
GMSH="${GMSH:-$(command -v gmsh || echo gmsh)}"
FAILS=0
fail() { echo "FAIL : $*"; FAILS=$((FAILS + 1)); }

mkdir -p "$WORK"
cd "$WORK" || exit 1
echo "date : $(date -Is) on $(hostname)"
echo "git  : $(git -C "$CBS3D_CODE" rev-parse --short HEAD) ($(git -C "$CBS3D_CODE" rev-parse --abbrev-ref HEAD))"

# ================================
# STEP 1: CASE
# ================================
"$GMSH" "$SRC/slab.geo" -3 -format msh2 -o slab.msh > gmsh.log 2>&1 || fail "gmsh"
python3 "$CBS3D_CODE/tools/gmsh_to_cbs.py" slab.msh --case-name slab \
    --out-dir serial --par "$SRC/slab.par" --matprop "$SRC/slab.matprop" \
    --overwrite > convert.log 2>&1 || fail "gmsh_to_cbs"
python3 "$CBS3D_CODE/tools/export_rank_local_cbs.py" serial/slab --parts 4 \
    --out-dir p4 --gmsh "$GMSH" --overwrite > export.log 2>&1 || fail "export 4 ranks"

# ================================
# STEP 2-3: SERIAL AND 4-RANK RUNS
# ================================
( cd serial && "$CBS3D_SERIAL_EXE" slab > solver.out 2>&1 ) || fail "serial run (see $WORK/serial/solver.out)"
mkdir -p mpi4
# ONE RANK PER PHYSICAL CORE. WITHOUT --hint/--cpu-bind THE 4 RANKS WERE
# BOUND TO 2 HYPERTHREADS AND RAN ~1000 TIMES SLOWER (BUSY-POLLING MPI).
# THE MPI RUN USES THE PRODUCTION SETTINGS OF blanket.sh (TEMPERATURE AFC,
# MONOTONE LTS). AFC LIMITS FLUID ELEMENTS ONLY, SO THE SOLID BAR MUST STILL
# EQUAL THE SERIAL GALERKIN FIELD TO ROUND-OFF.
( cd mpi4 && CBS3D_TEMPERATURE_AFC=1 CBS3D_LTS_MONOTONE=1 \
    srun --mpi="$CBS3D_SRUN_MPI" --ntasks=4 --cpus-per-task=1 \
    --hint=nomultithread --cpu-bind=cores --distribution=block:block \
    "$CBS3D_MPI_EXE" slab "$WORK/p4" \
    > solver.out 2>&1 ) || fail "4-rank run (see $WORK/mpi4/solver.out)"

# ================================
# STEP 4-5: COMPARE WITH THE EXACT SOLUTION
# ================================
for run in serial mpi4; do
    echo "---- $run ----"
    CMP=()
    [[ "$run" == mpi4 ]] && CMP=(--compare serial/output)     # MPI MUST EQUAL SERIAL
    python3 "$CBS3D_CODE/tools/check_nuclear_slab.py" "$run/output" "${CMP[@]}" || fail "$run temperature"
    Q=$(sed 's/\x1b\[[0-9;]*m//g' "$run/solver.out" | grep -oE "HEAT_INPUT_TOTAL_W = [-+0-9.eE]+" | tail -1 | awk '{print $3}')
    python3 - "$Q" <<'EOF' || fail "$run heat input"
import math, sys
q = float(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] else float("nan")
A, lam, L, W, qs = 1.6e6, 0.05, 0.10, 0.01, 1.0e5
exact = A * W * W * lam * (1.0 - math.exp(-L / lam)) + qs * W * W
rel = abs(q - exact) / exact
print("HEAT IN: printed %.6e W, exact %.6e W, rel. error %.2e" % (q, exact, rel))
sys.exit(0 if rel < 1.0e-3 else 1)
EOF
done

echo
if [[ "$FAILS" -eq 0 ]]; then
    echo "SELF-CHECK PASS"
else
    echo "SELF-CHECK FAIL ($FAILS problem(s)); FILES IN $WORK"
    exit 1
fi
