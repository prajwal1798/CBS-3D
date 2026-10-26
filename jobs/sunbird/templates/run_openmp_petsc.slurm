#!/bin/bash -l
# ================================
# FILE      : run_openmp_petsc.slurm
# PURPOSE   : ONE CBS3D PROCESS, OpenMP THREADS, PETSc AMG PRESSURE SOLVE.
#             RUNS A SERIAL CASE (.plt .bco .material .par .matprop).
#
# USAGE     : sbatch run_openmp_petsc.slurm CASE_DIR CASE_NAME
#             sbatch --time=24:00:00 run_openmp_petsc.slurm ...   (LONGER RUN)
#
# OUTPUT    : CASE_DIR/output/CASE_NAME/   (residuals.csv, .pvd, .vtu)
#             CASE_DIR/logs/               (solver log, run record)
#
# NOTES     : SET solver_opt = 3 IN THE .par ("Pressure Field Solver").
#             PETSc OPTIONS GO IN PETSC_OPTIONS: cbs3d_serial ACCEPTS ONLY
#             THE CASE NAME ON ITS COMMAND LINE.
# ================================
#SBATCH --job-name=cbs3d_omp
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --output=cbs3d_omp_%j.out
#SBATCH --error=cbs3d_omp_%j.err

set -euo pipefail

CASE_DIR="$(readlink -f "${1:?USAGE: sbatch run_openmp_petsc.slurm CASE_DIR CASE_NAME}")"
CASE="${2:?USAGE: sbatch run_openmp_petsc.slurm CASE_DIR CASE_NAME}"

source "$HOME/cbs3d_env.sh"

# ================================
# STEP 1: CHECK INPUTS
# ================================
for ext in plt bco material par matprop; do
    [[ -s "$CASE_DIR/$CASE.$ext" ]] || { echo "FATAL: missing $CASE_DIR/$CASE.$ext" >&2; exit 1; }
done
[[ -x "$CBS3D_SERIAL_EXE" ]] || { echo "FATAL: missing $CBS3D_SERIAL_EXE" >&2; exit 1; }

# ================================
# STEP 2: THREADS AND PETSc OPTIONS
# ================================
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export PETSC_OPTIONS="${PETSC_OPTIONS:--pc_type hypre -pc_hypre_type boomeramg}"

# ================================
# STEP 3: RUN RECORD
# ================================
mkdir -p "$CASE_DIR/logs"
RECORD="$CASE_DIR/logs/run_${SLURM_JOB_ID}.txt"
{
    echo "date        : $(date -Is)"
    echo "job         : $SLURM_JOB_ID on $(hostname)"
    echo "case        : $CASE_DIR/$CASE"
    echo "threads     : $OMP_NUM_THREADS"
    echo "PETSC_OPTIONS: $PETSC_OPTIONS"
    echo "exe         : $CBS3D_SERIAL_EXE"
    echo "git commit  : $(git -C "$CBS3D_CODE" rev-parse HEAD 2>/dev/null || echo unknown)"
    sha256sum "$CBS3D_SERIAL_EXE" "$CASE_DIR/$CASE".{plt,bco,material,par,matprop}
} > "$RECORD"

# ================================
# STEP 4: RUN (FROM THE CASE DIRECTORY; THE SOLVER WRITES output/<case>)
# ================================
cd "$CASE_DIR"
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    stdbuf -oL -eL "$CBS3D_SERIAL_EXE" "$CASE" \
    > "logs/solver_${SLURM_JOB_ID}.out" 2> "logs/solver_${SLURM_JOB_ID}.err"

sacct -j "$SLURM_JOB_ID" --format=JobID,Elapsed,State,ExitCode,MaxRSS >> "$RECORD" || true
echo "DONE: $CASE_DIR/output/$CASE"
