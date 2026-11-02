#!/bin/bash -l
# G8b 2-vs-4 MPI partition consistency, actual U-channel hot step-50000.
# Creates isolated 4-rank METIS partition and remaps existing 2-rank fields
# by GLOBAL ID. No interpolation, original inputs or production code edits.
#SBATCH --job-name=cbs3d_g8b
#SBATCH --account=scw7008
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=3G
#SBATCH --time=00:20:00
#SBATCH --output=cbs3d_g8b_%j.out
#SBATCH --error=cbs3d_g8b_%j.err

set -euo pipefail
source "$HOME/cbs3d_env.sh"
export OMP_NUM_THREADS=1
BASE="$CBS3D_SCRATCH/cases/uchannel_afc_audit"
ORIGINAL="$BASE/run_uchannel_29706"
G6ROOT="$CBS3D_SCRATCH/verification_g6/29763"
REF="$G6ROOT/run_uchannel_29763/output/uchannel"
MASTER="$G6ROOT/flow_one_step.par"
PAR="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.par"
MAT="$ORIGINAL/partitions/rank_0000/uchannel_rank_0000.matprop"
WORK="$CBS3D_SCRATCH/verification_g8b/${SLURM_JOB_ID:-manual}"
PART4="$WORK/partitions_4"
REMAP="$WORK/restart_vtu_step_50000_four_rank"
mkdir -p "$WORK"

check_hash() {
    local actual
    actual="$(sha256sum "$2" | awk '{print $1}')"
    echo "SHA256 $2 $actual"
    [[ "$actual" == "$1" ]] || { echo "G8b FATAL original hash mismatch: $2"; exit 1; }
}
for f in "$MASTER" "$PAR" "$MAT" \
    "$REF/uchannel_step_00050000_rank_0000.vtu" \
    "$REF/uchannel_step_00050000_rank_0001.vtu" \
    "$REF/uchannel_step_00050001_rank_0000.vtu" \
    "$REF/uchannel_step_00050001_rank_0001.vtu"; do
    [[ -s "$f" ]] || { echo "G8b FATAL missing source: $f"; exit 1; }
done
check_hash "96ca7134f6affc2cd17cdfb6b14ba70bbda9f243d7041ba196d70a00bc1c16fa" "$CBS3D_MPI_EXE"
check_hash "706b64a8e641a53294c393c23dba3c086cd104fd2961541af7cfe3d01839e6ac" "$PAR"
check_hash "928bc80463893a6d80c19eb9b5df27e5553ccd8b9821473935acf7e12641a3fb" "$MAT"
cmp -s "$MASTER" "$G6ROOT/run_uchannel_29763/partitions/rank_0000/uchannel_rank_0000.par" \
    || { echo "G8b FATAL reference one-step .par changed"; exit 1; }

# Locate the original UNPARTITIONED CBS .plt/.bco/.material triple.
# Avoid accidentally exporting the old two-rank mesh (which has fewer cells).
SERIAL="$(python3 - "$BASE" "$CBS3D_SCRATCH/cases" "$CBS3D_CODE" <<'PY'
import os
import sys
from pathlib import Path
override = os.environ.get("CBS3D_G8B_SERIAL_BASE","").strip()
if override:
    bases = [Path(override)]
else:
    bases = []
    roots = [Path(x) for x in sys.argv[1:]]
    for root in roots:
        if not root.is_dir():
            continue
        for folder,dirs,files in os.walk(root):
            depth = len(Path(folder).relative_to(root).parts)
            dirs[:] = [
                d for d in dirs if depth < 3 and
                not (d.startswith(("run_","partitions_","verification_","output","build")))
            ]
            for name in files:
                if name.endswith(".plt"):
                    bases.append(Path(folder)/name[:-4])
found = []
for base in set(bases):
    if not all(Path(str(base)+ext).is_file() for ext in (".plt",".bco",".material")):
        continue
    try:
        with open(str(base)+".plt",encoding="utf-8") as f:
            header=list(map(int,f.readline().split()))
    except (OSError,ValueError):
        continue
    if header[:2] == [45567,9613] and header[2] == 8860:
        found.append(base)
if len(found)!=1:
    sys.stderr.write(f"G8b mesh discovery: expected exactly ONE full U-channel serial base, got {len(found)}:\n")
    for candidate in sorted(found):
        sys.stderr.write(f"  {candidate}\n")
    sys.stderr.write("Set CBS3D_G8B_SERIAL_BASE=/path/to/serial/uchannel (without extension) if needed.\n")
    raise SystemExit(1)
print(found[0])
PY
)"
echo "G8b original serial mesh: $SERIAL"

# Partition EXACT same global mesh anew with four independent ranks.
python3 "$CBS3D_CODE/tools/export_rank_local_cbs.py" "$SERIAL" \
    --parts 4 \
    --out-dir "$PART4" \
    --case-name uchannel \
    --par "$MASTER" \
    --matprop "$MAT"

# Validate partition inventory. No mesh rescaling, no new nodes or elements.
python3 - "$PART4/partition_manifest.json" <<'PY'
import json,sys
m=json.load(open(sys.argv[1],encoding="utf-8"))
assert m["parts"]==4 and m["global"]=={
    "tetrahedra":45567,"nodes":9613,"physical_boundary_triangles":8860
}, m["global"]
assert all(m["checks"].values()),m["checks"]
print("G8b 4-rank partition topology check PASS:",m["global"])
PY

# Generate rank-local step-50000 restart VTU files using the *original*
# 2-rank OWNED global nodal fields; no interpolation or serial recompute.
python3 "$CBS3D_CODE/tools/verify_uchannel_partition_parity_g8b.py" make-restart \
    --reference "$REF" --partitions "$PART4" --destination "$REMAP"

# Isolate optional environment switches to match Job29763 exactly.
for name in \
    CBS3D_AVERAGE_START_ITERATION CBS3D_AVERAGE_START_TIME \
    CBS3D_CHECKPOINT_ROOT CBS3D_CONTINUATION_ITERATIONS \
    CBS3D_RESTART_VELOCITY_ARRAY CBS3D_SA_WALL_TREATMENT \
    CBS3D_CHT_WALL_TREATMENT \
    CBS3D_STEADY_ENERGY_MAXIT CBS3D_STEADY_ENERGY_RELAX \
    CBS3D_STEADY_ENERGY_TOL CBS3D_STEADY_ENERGY_RESIDUAL_TOL
do
    unset "$name"
done
export CBS3D_TEMPERATURE_AFC=0
export CBS3D_LTS_MONOTONE=1
export CBS3D_LTS_REFRESH_MARGIN=0.8
export CBS3D_STEADY_ENERGY=0
export CBS3D_CHECKPOINT_EVERY=0
export CBS3D_RESTART=1
export CBS3D_RESTART_FORMAT=legacy_vtu
export CBS3D_RESTART_ROOT="$REMAP"
export CBS3D_RESTART_ITERATION=50000
export CBS3D_RESTART_AS_INITIAL_CONDITION=0
export CBS3D_PAR="$MASTER"
export CBS3D_MATPROP="$MAT"

echo "DATE: $(date -Is)"
echo "G8b source: $(git -C "$CBS3D_CODE" rev-parse HEAD)"
echo "Work: $WORK"
echo "2-rank original replay: $REF"
echo "4-rank remapped initial state: $REMAP"
export SLURM_SUBMIT_DIR="$WORK"
bash "$CBS3D_CODE/jobs/sunbird/templates/run_mpi.slurm" "$PART4" uchannel

CANDIDATE="$WORK/run_uchannel_${SLURM_JOB_ID:-manual}/output/uchannel"
python3 "$CBS3D_CODE/tools/verify_uchannel_partition_parity_g8b.py" compare \
    --reference "$REF" \
    --candidate "$CANDIDATE"
echo "G8B PARTITION-COUNT PARITY PASS"
