#!/bin/bash -l
# ================================
# FILE      : run_mpi_anycores.slurm
# PURPOSE   : SAME RUN AS run_mpi.slurm, BUT SLURM MAY PLACE THE N RANKS ON
#             ANY FREE CORES OF ANY NODES (E.G. 20 ON ONE NODE, 30 ON ANOTHER,
#             128 ON A THIRD). NO --nodes, NO --exclusive, MEMORY PER CORE.
#
# USAGE     : sbatch --ntasks=N run_mpi_anycores.slurm PARTITION_ROOT CASE_NAME
#             (PARTITION_ROOT MUST HOLD EXACTLY N rank_XXXX FOLDERS)
#             EASIER: bash submit_mpi.sh CASE_DIR CASE N any
#
# NOTES     : USE WHEN THE CLUSTER IS BUSY. WHEN WHOLE NODES ARE FREE,
#             run_mpi.slurm (EXCLUSIVE NODES) IS FASTER: NO SHARED MEMORY
#             BANDWIDTH WITH OTHER JOBS, FEWER NODES TO TALK BETWEEN.
#             THE SLOWEST RANK SETS THE PACE OF EVERY ITERATION.
# ================================
#SBATCH --job-name=cbs3d_mpi_any
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --output=cbs3d_mpi_%j.out
#SBATCH --error=cbs3d_mpi_%j.err

set -euo pipefail
source "$HOME/cbs3d_env.sh"
exec bash "$CBS3D_CODE/jobs/sunbird/templates/run_mpi.slurm" "$@"
