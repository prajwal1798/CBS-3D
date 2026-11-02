#!/bin/bash -l
# ================================
# FILE      : verify_two_layer_g4_no_afc.slurm
# PURPOSE   : TWO-MATERIAL CONFORMAL FLUID/SOLID CONDUCTION PATCH TEST.
#             EXACT PIECEWISE-LINEAR TEMPERATURE IN y:
#               0. CLEAR PRODUCTION CONTROLS INHERITED FROM THE SUBMITTING
#                  SHELL (RESTART, AVERAGES, CHECKPOINTS, STEADY ENERGY)
#               1. MESH cht2.geo, CONVERT TO A SERIAL CBS CASE
#               2. SERIAL RUN    (CBS3D_SERIAL_EXE)
#               3. 4-RANK MPI RUN (CBS3D_MPI_EXE ON 4 PARTITIONS)
#                  (WITH CBS3D_TEMPERATURE_AFC=0, LIKE SERIAL)
#               4. REQUIRE STEADY STOP, THEN BOTH AGAINST EXACT TEMPERATURE (tools/check_two_layer_conduction.py),
#                  AND THE MPI FIELD AGAINST THE SERIAL FIELD, NODE BY NODE
#               5. THE PRINTED HEAT INPUT AGAINST q''*W^2=0.016 W
#             LAST LINE: G4-CONFORMAL-NO-AFC PASS OR G4-CONFORMAL-NO-AFC FAIL.
#
# USAGE     : sbatch verify_two_layer_g4_no_afc.slurm          (selfcheck.sh SUBMITS IT)
#
# OUTPUT    : $CBS3D_SCRATCH/verification_g4/<JOBID>/ AND cbs3d_g4_<JOBID>.out
# ================================
#SBATCH --job-name=cbs3d_g4
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=cbs3d_g4_%j.out
#SBATCH --error=cbs3d_g4_%j.err

set -uo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1

# ================================
# HERMETIC ENVIRONMENT
# ================================
# blanket.sh run EXPORTS RESTART/AVERAGE/CHECKPOINT/STEADY-ENERGY CONTROLS,
# AND sbatch PASSES THE WHOLE SUBMITTING ENVIRONMENT TO THIS JOB. JOB 29536
# INHERITED THE L1T RESTART AND TRIED TO OPEN
#   run_blanket_29513/output/blanket/slab_step_00100000_rank_*.vtu
# SO THE 4-RANK STEP DIED BEFORE WRITING mpi4/output. THE SELF-CHECK MUST
# DEPEND ONLY ON THIS FILE AND THE CASE, NEVER ON THE SUBMITTING SHELL.
for variable in \
    CBS3D_AFC_LOG_EVERY CBS3D_AFC_MIN_TEMPERATURE CBS3D_AFC_SOLID \
    CBS3D_AVERAGE_START_ITERATION CBS3D_AVERAGE_START_TIME \
    CBS3D_CHECKPOINT_EVERY CBS3D_CHECKPOINT_ROOT \
    CBS3D_CHT_WALL_TREATMENT CBS3D_CONTINUATION_ITERATIONS \
    CBS3D_LTS_REFRESH_MARGIN CBS3D_MATPROP CBS3D_PAR \
    CBS3D_RESTART CBS3D_RESTART_AS_INITIAL_CONDITION \
    CBS3D_RESTART_FORMAT CBS3D_RESTART_ITERATION CBS3D_RESTART_ROOT \
    CBS3D_RESTART_VELOCITY_ARRAY CBS3D_SA_WALL_TREATMENT \
    CBS3D_STEADY_ENERGY CBS3D_STEADY_ENERGY_MAXIT \
    CBS3D_STEADY_ENERGY_RELAX CBS3D_STEADY_ENERGY_TOL \
    CBS3D_STEADY_ENERGY_RESIDUAL_TOL \
    CBS3D_STEADY_ENERGY_ANDERSON_M \
    CBS3D_VERBOSE CBS3D_WALL_B CBS3D_WALL_KAPPA
do
    unset "$variable"
done

# Both solvers run without AFC and without the experimental implicit energy path.
export CBS3D_TEMPERATURE_AFC=0
export CBS3D_LTS_MONOTONE=0

SRC="$CBS3D_CODE/cases/two_layer_conduction"
WORK="$CBS3D_SCRATCH/verification_g4/${SLURM_JOB_ID:-local}"
GMSH="${GMSH:-$(command -v gmsh || echo gmsh)}"
FAILS=0
fail() { echo "FAIL : $*"; FAILS=$((FAILS + 1)); }

mkdir -p "$WORK"
cd "$WORK" || exit 1
echo "date : $(date -Is) on $(hostname)"
echo "git  : $(git -C "$CBS3D_CODE" rev-parse --short HEAD) ($(git -C "$CBS3D_CODE" rev-parse --abbrev-ref HEAD))"
echo "env  : AFC=${CBS3D_TEMPERATURE_AFC} LTS=${CBS3D_LTS_MONOTONE} restart=${CBS3D_RESTART:-0} steady=${CBS3D_STEADY_ENERGY:-0} average=${CBS3D_AVERAGE_START_TIME:-off}"
sha256sum "$CBS3D_SERIAL_EXE" "$CBS3D_MPI_EXE"

# ================================
# STEP 1: CASE
# ================================
"$GMSH" "$SRC/cht2.geo" -3 -format msh2 -o cht2.msh > gmsh.log 2>&1 || fail "gmsh"
python3 "$CBS3D_CODE/tools/gmsh_to_cbs.py" cht2.msh --case-name cht2 \
    --out-dir serial --par "$SRC/cht2_convergence.par" --matprop "$SRC/cht2.matprop" \
    --overwrite > convert.log 2>&1 || fail "gmsh_to_cbs"
python3 "$CBS3D_CODE/tools/export_rank_local_cbs.py" serial/cht2 --parts 4 \
    --out-dir p4 --gmsh "$GMSH" --overwrite > export.log 2>&1 || fail "export 4 ranks"

# ================================
# STEP 2-3: SERIAL AND 4-RANK RUNS
# ================================
( cd serial && "$CBS3D_SERIAL_EXE" cht2 > solver.out 2>&1 ) || fail "serial run (see $WORK/serial/solver.out)"
mkdir -p mpi4
# One rank per physical core. Both runs use AFC=0.
# Partition independence is assessed using the same exact cht2 problem.
( cd mpi4 && CBS3D_TEMPERATURE_AFC=0 CBS3D_LTS_MONOTONE=0 \
    srun --mpi="$CBS3D_SRUN_MPI" --ntasks=4 --cpus-per-task=1 \
    --hint=nomultithread --cpu-bind=cores --distribution=block:block \
    "$CBS3D_MPI_EXE" cht2 "$WORK/p4" \
    > solver.out 2>&1 ) || fail "4-rank run (see $WORK/mpi4/solver.out)"

# ================================
# STEP 4-5: COMPARE WITH THE EXACT SOLUTION
# ================================
for run in serial mpi4; do
    echo "---- $run ----"
    # The comparison is against a STATIONARY analytic solution; the numerical
    # run is invalid until the thermal steady-state convergence gate is met.
    # Catch an iteration-limit stop even if an unconverged field happens to
    # fall within the benchmark's spatial error tolerance.
    if ! grep -qi "steady-state convergence reached" "$run/solver.out"; then
        fail "$run terminated without steady-state convergence"
        grep -i -E "stop reason|maximum iteration count" "$run/solver.out" | tail -n 3 || true
    fi
    CMP=()
    [[ "$run" == mpi4 ]] && CMP=(--compare serial/output)     # MPI MUST EQUAL SERIAL
    python3 "$CBS3D_CODE/tools/check_two_layer_conduction.py" "$run/output" "${CMP[@]}" || fail "$run temperature"
    Q=$(sed 's/\x1b\[[0-9;]*m//g' "$run/solver.out" | grep -oE "HEAT_INPUT_TOTAL_W = [-+0-9.eE]+" | tail -1 | awk '{print $3}')
    python3 - "$Q" <<'EOF' || fail "$run heat input"
import sys
q = float(sys.argv[1]) if len(sys.argv)>1 and sys.argv[1] else float("nan")
expected = 1000.0 * 0.004 * 0.004
rel = abs(q - expected) / expected
print("HEAT IN: printed %.9e W, exact %.9e W, rel. error %.2e" % (q, expected, rel))
sys.exit(0 if rel < 1.0e-3 else 1)
EOF
done

echo
if [[ "$FAILS" -eq 0 ]]; then
    echo "G4-CONFORMAL-NO-AFC PASS"
else
    echo "G4-CONFORMAL-NO-AFC FAIL ($FAILS problem(s)); FILES IN $WORK"
    exit 1
fi
