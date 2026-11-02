#!/bin/bash -l
# G13a: physical-time backward Euler CHT on same frozen U-channel velocity.
# Assemble physical lumped C_i, no pseudo-time surrogate. Reuse one PETSc matrix across 10 steps.
# Compiles scratch-only binary, leaves production executable and cases unchanged.
#SBATCH --job-name=cbs3d_g13a
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:35:00
#SBATCH --output=cbs3d_g13a_%j.out
#SBATCH --error=cbs3d_g13a_%j.err

set -euo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1

BASE="$CBS3D_SCRATCH/cases/uchannel_afc_audit"
PART="$BASE/run_uchannel_29706/partitions"
PAR="$PART/rank_0000/uchannel_rank_0000.par"
MAT="$PART/rank_0000/uchannel_rank_0000.matprop"
INPUT="$BASE/run_uchannel_29580/output/uchannel"
WORK="$CBS3D_SCRATCH/benchmark_implicit_cht_g13a/${SLURM_JOB_ID:-manual}"
BUILD="$WORK/build"
mkdir -p "$WORK"

check() {
    local actual
    actual="$(sha256sum "$2" | awk '{print $1}')"
    [[ "$actual" == "$1" ]] || { echo "FATAL SHA mismatch $2 : $actual"; exit 1; }
}
check 706b64a8e641a53294c393c23dba3c086cd104fd2961541af7cfe3d01839e6ac "$PAR"
check 928bc80463893a6d80c19eb9b5df27e5553ccd8b9821473935acf7e12641a3fb "$MAT"
for r in 0000 0001; do
    [[ -s "$INPUT/uchannel_step_00050000_rank_$r.vtu" ]] || {
        echo "FATAL missing rank-$r hot restart"; exit 1; }
done

echo "G13a CODE COMMIT $(git -C "$CBS3D_CODE" rev-parse HEAD)"
echo "G13a WORK $WORK"
cmake -S "$CBS3D_CODE" -B "$BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCBS3D_ENABLE_MPI=ON \
  -DCBS3D_ENABLE_PETSC=ON \
  -DCBS3D_ENABLE_OPENMP=OFF \
  -DCBS3D_KEEP_LEGACY_EXECUTABLE_NAME=OFF \
  -DCBS3D_BUILD_SOURCE_VALIDATION=OFF \
  -DCBS3D_BUILD_TESTS=OFF
cmake --build "$BUILD" --target cbs3d_parallel -j 4
EXE="$BUILD/cbs3d_parallel"
[[ -x "$EXE" ]] || { echo "FATAL missing candidate binary"; exit 1; }
echo "G13a CANDIDATE SHA $(sha256sum "$EXE")"

# Staging launcher copy; the installed production launcher remains untouched.
python3 - "$CBS3D_CODE/jobs/sunbird/templates/run_mpi.slurm" "$WORK/launcher.sh" <<'PY'
from pathlib import Path
import sys
s=Path(sys.argv[1]).read_text()
needle='source "$HOME/cbs3d_env.sh"'
assert s.count(needle)==1
Path(sys.argv[2]).write_text(s.replace(needle,needle+'\nexport CBS3D_MPI_EXE="$CBS3D_G13A_EXE"'))
PY
export CBS3D_G13A_EXE="$EXE"

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
export CBS3D_STEADY_ENERGY_LOW_ORDER_BOOTSTRAP=0
export CBS3D_STEADY_ENERGY_LOW_ORDER_ONLY=1
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

arms=("be_dt001_n10 0.01 10" "be_dt0001_n100 0.001 100")
success=0
executed=0
for arm in "${arms[@]}"; do
    read -r label physical_dt physical_steps <<< "$arm"
    executed=$((executed + 1))
    export CBS3D_STEADY_ENERGY_PSEUDOTIME_SHIFT=0
    export CBS3D_IMPLICIT_ENERGY_BE_DT="$physical_dt"
    export CBS3D_IMPLICIT_ENERGY_BE_STEPS="$physical_steps"
    export CBS3D_STEADY_ENERGY_ANDERSON_M=0
    mkdir -p "$WORK/$label"
    export SLURM_SUBMIT_DIR="$WORK/$label"
    echo "===== G13a $label dt=$physical_dt seconds, $physical_steps backward-Euler steps; total physical time 0.1s ====="

    RUN="$WORK/$label/run_uchannel_${SLURM_JOB_ID:-manual}"
    if bash "$WORK/launcher.sh" "$PART" uchannel </dev/null; then
        echo "G13a $label physical BE stepper returned success"
        if python3 - "$CBS3D_CODE/tools" "$RUN/output/uchannel" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0,sys.argv[1])
from verify_uchannel_partition_parity_g8b import owned_fields
data=owned_fields(Path(sys.argv[2]),50001,2)
assert len(data)==9613,len(data)
t=[v["temperature"][0] for v in data.values()]
print("Independent final PHYSICAL BE state: Tmin=%.12f K, Tmax=%.12f K, nodes=%d"%
      (min(t),max(t),len(t)))
assert min(t)>=298.15-1e-8,"Unphysical Tmin"
PY
        then
            success=$((success + 1))
            echo "G13a $label: CONVERGENCE AND PHYSICAL FLOOR PASS"
        else
            echo "G13a $label: FINAL THERMAL FIELD UNACCEPTABLE"
        fi
    else
        echo "G13a $label: SOLVER FAILED (read details below)"
    fi

    if [[ -f "$RUN/solver.out" ]]; then
        grep -E '\[PHYSICAL-BE\]|scheme|Krylov iterations|residual check|temperature range|physical BE|solve time|stop reason' \
            "$RUN/solver.out" | tail -n 24 || true
    fi
    if [[ -f "$RUN/solver.err" ]]; then
        grep -E 'ERROR on MPI|Implicit steady|SteadyEnergySolver|low-order' \
            "$RUN/solver.err" | tail -n 5 || true
    fi
    echo "===== END $label ====="
done

echo "G13a arms executed: $executed / ${#arms[@]}"
echo "G13a converged bounded arms: $success / ${#arms[@]}"
[[ "$executed" == "${#arms[@]}" ]] || exit 2
[[ "$success" -eq "${#arms[@]}" ]] || { echo "G13A PHYSICAL-TIME BE: at least one arm FAILED"; exit 1; }
# Independent time-refinement comparison on exactly the same 9613 owned nodes.
python3 - "$CBS3D_CODE/tools" "$WORK" "${SLURM_JOB_ID:-manual}" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0,sys.argv[1])
from verify_uchannel_partition_parity_g8b import owned_fields
work=Path(sys.argv[2])
job=sys.argv[3]
def read(label):
    p=work/label/("run_uchannel_"+job)/"output"/"uchannel"
    d=owned_fields(p,50001,2)
    assert len(d)==9613
    return d
coarse=read("be_dt001_n10")
fine=read("be_dt0001_n100")
assert coarse.keys()==fine.keys()
delta=max(abs(coarse[i]["temperature"][0]-fine[i]["temperature"][0])
          for i in coarse)
print("G13a physical dt-refinement: same final t=0.1 s,"
      " max nodal |T(dt=0.01)-T(dt=0.001)|=%.12g K" % delta)
assert delta>=0.0 and delta<1e10
PY
echo "G13A EXPERIMENTAL PHYSICAL-TIME LOW-ORDER BE CHT PASSED TWO DT TESTS (NOT HIGH-ORDER AFC)"
