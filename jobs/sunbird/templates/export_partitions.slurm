#!/bin/bash -l
# ================================
# FILE      : export_partitions.slurm
# PURPOSE   : SPLIT A VALIDATED SERIAL CASE INTO N RANK-LOCAL CASES ON A
#             COMPUTE NODE (LARGE MESHES NEED SEVERAL GB; NOT ON THE LOGIN NODE).
#
# USAGE     : sbatch export_partitions.slurm SERIAL_BASE NPARTS OUT_DIR
#   SERIAL_BASE  e.g. /scratch/$USER/CBS3D/cases/fw/serial/fw  (NO EXTENSION)
#
# NEEDS     : A gmsh EXECUTABLE (MODULE OR $HOME/.local/bin/gmsh FROM pip).
# ================================
#SBATCH --job-name=cbs3d_part
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=cbs3d_part_%j.out
#SBATCH --error=cbs3d_part_%j.err

set -euo pipefail

BASE="${1:?USAGE: sbatch export_partitions.slurm SERIAL_BASE NPARTS OUT_DIR}"
NPARTS="${2:?NPARTS}"
OUT="${3:?OUT_DIR}"

source "$HOME/cbs3d_env.sh"
GMSH="${GMSH:-$(command -v gmsh || echo "$HOME/.local/bin/gmsh")}"

python3 "$CBS3D_CODE/tools/mesh_summary.py" "$BASE" --quiet
python3 "$CBS3D_CODE/tools/export_rank_local_cbs.py" "$BASE" \
    --parts "$NPARTS" --out-dir "$OUT" --gmsh "$GMSH"
python3 "$CBS3D_CODE/tools/dd4d_check_partition.py" "$OUT" "$(basename "$BASE")" "$NPARTS"
