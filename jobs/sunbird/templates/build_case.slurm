#!/bin/bash -l
# ================================
# FILE      : build_case.slurm
# PURPOSE   : tools/cbs_case.sh ON A COMPUTE NODE: GMSH MESH, CBS CONVERSION,
#             MESH CHECK AND (OPTIONALLY) N RANK-LOCAL PARTITIONS.
#             USE IT FOR MESHES ABOVE ~0.5 M TETS (MINUTES AND SEVERAL GB).
#
# USAGE     : sbatch build_case.slurm GEO PAR MATPROP OUT_DIR [cbs_case.sh OPTIONS]
#   e.g.      sbatch $T/build_case.slurm \
#                 $CBS3D_CODE/cases/lid_driven_cavity/ldc3d.geo \
#                 $CBS3D_CODE/cases/lid_driven_cavity/ldc3d.par \
#                 $CBS3D_CODE/cases/lid_driven_cavity/ldc3d.matprop \
#                 $CBS3D_SCRATCH/cases/ldc3d --parts 128 --overwrite
#
# OUTPUT    : OUT_DIR/{mesh,serial,partitions_N}  (SEE tools/cbs_case.sh)
#             cbs3d_case_<JOBID>.out IN THE SUBMIT DIRECTORY
# ================================
#SBATCH --job-name=cbs3d_case
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=cbs3d_case_%j.out
#SBATCH --error=cbs3d_case_%j.err

set -euo pipefail
[[ "$#" -ge 4 ]] || { echo "USAGE: sbatch build_case.slurm GEO PAR MATPROP OUT_DIR [OPTIONS]" >&2; exit 1; }

source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

echo "date : $(date -Is) on $(hostname)"
echo "git  : $(git -C "$CBS3D_CODE" rev-parse --short HEAD)"
bash "$CBS3D_CODE/tools/cbs_case.sh" "$@" \
    --python "$(command -v python3)" --threads "$SLURM_CPUS_PER_TASK"
echo "elapsed: ${SECONDS} s"
sacct -j "$SLURM_JOB_ID" --format=JobID,Elapsed,State,MaxRSS || true
