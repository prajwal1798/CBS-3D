#!/bin/bash -l
# G12b: read-only controlled sensitivity of existing implicit AFC solver.
# Reuse *identical* isolated corrected binary from Job 29833; DO NOT rebuild.
# Vary only AFC preconditioner pseudo-time shift and Anderson depth.
# Failures are recorded per arm without aborting other arms.
#SBATCH --job-name=cbs3d_g12b
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=3G
#SBATCH --time=00:20:00
#SBATCH --output=cbs3d_g12b_%j.out
#SBATCH --error=cbs3d_g12b_%j.err

set -euo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1
BASE="$CBS3D_SCRATCH/cases/uchannel_afc_audit"
PART="$BASE/run_uchannel_29706/partitions"
PAR="$PART/rank_0000/uchannel_rank_0000.par"
MAT="$PART/rank_0000/uchannel_rank_0000.matprop"
INPUT="$BASE/run_uchannel_29580/output/uchannel"
EXE="$CBS3D_SCRATCH/benchmark_implicit_cht_g12/29833/build/cbs3d_parallel"
WORK="$CBS3D_SCRATCH/benchmark_implicit_cht_g12b/${SLURM_JOB_ID:-manual}"
mkdir -p "$WORK"

[[ -x "$EXE" ]] || { echo "G12b FATAL corrected Job29833 binary missing"; exit 1; }
for f in "$PAR" "$MAT" \
    "$INPUT/uchannel_step_00050000_rank_0000.vtu" \
    "$INPUT/uchannel_step_00050000_rank_0001.vtu"; do
    [[ -s "$f" ]] || { echo "G12b FATAL missing input: $f"; exit 1; }
done
verify() {
    local actual
    actual="$(sha256sum "$2" | awk '{print $1}')"
    echo "SHA256 $2 $actual"
    [[ "$actual" == "$1" ]] || { echo "G12b FATAL checksum mismatch: $2"; exit 1; }
}
verify "6a565fb90209b0f8aa15a5b3c8ab4919cfff2f8e79d10d9a70321ae5c7de0f0d" "$EXE"
verify "706b64a8e641a53294c393c23dba3c086cd104fd2961541af7cfe3d01839e6ac" "$PAR"
verify "928bc80463893a6d80c19eb9b5df27e5553ccd8b9821473935acf7e12641a3fb" "$MAT"

python3 - "$CBS3D_CODE/jobs/sunbird/templates/run_mpi.slurm" "$WORK/launcher.sh" <<'PY'
from pathlib import Path
import sys
s=Path(sys.argv[1]).read_text()
hook='source "$HOME/cbs3d_env.sh"'
assert s.count(hook)==1
Path(sys.argv[2]).write_text(s.replace(hook,hook+'\nexport CBS3D_MPI_EXE="$CBS3D_G12B_EXE"'))
PY
export CBS3D_G12B_EXE="$EXE"

# Clean unrelated overrides: same mesh, flow restart, and thermal source.
for name in \
    CBS3D_AVERAGE_START_ITERATION CBS3D_AVERAGE_START_TIME \
    CBS3D_CHECKPOINT_ROOT CBS3D_CONTINUATION_ITERATIONS \
    CBS3D_RESTART_VELOCITY_ARRAY CBS3D_SA_WALL_TREATMENT \
    CBS3D_CHT_WALL_TREATMENT CBS3D_AFC_SOLID \
    CBS3D_AFC_MIN_TEMPERATURE \
    CBS3D_STEADY_ENERGY_ARMIJO_C1
do unset "$name"; done
export CBS3D_TEMPERATURE_AFC=1
export CBS3D_STEADY_ENERGY=1
export CBS3D_RESTART=1
export CBS3D_RESTART_FORMAT=legacy_vtu
export CBS3D_RESTART_ROOT="$INPUT"
export CBS3D_RESTART_ITERATION=50000
export CBS3D_RESTART_AS_INITIAL_CONDITION=0
export CBS3D_CHECKPOINT_EVERY=0
export CBS3D_PAR="$PAR"
export CBS3D_MATPROP="$MAT"
export CBS3D_LTS_MONOTONE=1
export CBS3D_LTS_REFRESH_MARGIN=0.8
export CBS3D_STEADY_ENERGY_MAXIT=150
export CBS3D_STEADY_ENERGY_TOL=1e-6
export CBS3D_STEADY_ENERGY_RESIDUAL_TOL=1e-6
export CBS3D_STEADY_ENERGY_RELAX=0.5
export CBS3D_AFC_LOG_EVERY=10000
export PETSC_OPTIONS="-energy_ksp_type gmres -energy_pc_type asm -energy_sub_pc_type ilu -energy_sub_pc_factor_levels 1 -energy_ksp_rtol 1e-10 -energy_ksp_max_it 4000"

echo "G12b nonlinear preconditioner parameter sensitivity; no solver source changes"
echo "CODE $(git -C "$CBS3D_CODE" rev-parse HEAD)"
echo "WORK $WORK"

# Label, pseudo-time coefficient gamma, Anderson depth:
# gamma=0 is exact steady low-order preconditioner; gamma=1 corresponds
# to C/dt and was too close to explicit relaxation in Job29833.
arm_table=$(cat <<'EOF'
shift1_anderson3 1 3
shift001 0.01 0
shift0001 0.0001 0
shift0 0 0
shift001_anderson3 0.01 3
EOF
)
success=0
while read -r label gamma anderson; do
    [[ -n "$label" ]] || continue
    export CBS3D_STEADY_ENERGY_PSEUDOTIME_SHIFT="$gamma"
    export CBS3D_STEADY_ENERGY_ANDERSON_M="$anderson"
    mkdir -p "$WORK/$label"
    export SLURM_SUBMIT_DIR="$WORK/$label"
    echo "===== G12b ARM $label gamma=$gamma Anderson=$anderson ====="
    if bash "$WORK/launcher.sh" "$PART" uchannel; then
        success=$((success+1))
        echo "G12b ARM $label CONVERGED"
    else
        echo "G12b ARM $label FAILED (nonconvergence or solver failure)"
    fi
    RUN="$WORK/$label/run_uchannel_${SLURM_JOB_ID:-manual}"
    if [[ -f "$RUN/solver.out" ]]; then
        grep -E '\[AFC\] iteration (1:|5:|10:|20:|50:|100:|150:)|Steady energy \(implicit|Krylov iterations|explicit residual|temperature range|solve time|stop reason' "$RUN/solver.out" | tail -n 20 || true
    fi
    if [[ -f "$RUN/solver.err" ]]; then
        grep -E 'ERROR on MPI|Implicit steady|non-finite|SteadyEnergySolver' "$RUN/solver.err" | tail -n 4 || true
    fi
    echo "===== END $label ====="
done <<< "$arm_table"

echo "G12b converged arms: $success / 5"
echo "G12B SENSITIVITY SWEEP COMPLETED (nonconverged arms remain failures)"
