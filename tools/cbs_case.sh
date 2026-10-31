#!/usr/bin/env bash
# ================================
# FILE      : cbs_case.sh
# PURPOSE   : ONE COMMAND FROM A GMSH .geo FILE TO A READY CBS3D CASE.
#               1. GMSH CLI       .geo  -> .msh (MSH 2.2)
#                  (OR GIVE A READY .msh INSTEAD OF THE .geo: IT IS COPIED
#                   AND NOT RE-MESHED; --clscale IS THEN IGNORED)
#               2. gmsh_to_cbs.py .msh  -> .plt .bco .material .par .matprop
#               3. mesh_summary.py      -> RESULT: PASS / FAIL
#               4. (--parts N)          -> N RANK-LOCAL CASES FOR MPI
#                  export_rank_local_cbs.py + dd4d_check_partition.py
#
# USAGE     : bash tools/cbs_case.sh GEO PAR MATPROP OUT_DIR [OPTIONS]
#
#   GEO       GMSH GEOMETRY WITH Physical Surface(...)/Physical Volume(...)
#             OR A MESH (.msh, ANY GMSH FORMAT) MADE FROM SUCH A GEOMETRY
#   PAR       SOLVER CONTROLS (COPIED UNCHANGED)
#   MATPROP   MATERIAL TABLE (COPIED UNCHANGED)
#   OUT_DIR   NEW CASE FOLDER. LAYOUT:
#               OUT_DIR/mesh/NAME.msh          GMSH MESH + gmsh.log
#               OUT_DIR/serial/NAME.*          SERIAL CASE (RUN FROM HERE)
#               OUT_DIR/partitions_N/rank_XXXX MPI CASE    (WITH --parts N)
#
# OPTIONS   : --name NAME     CASE NAME           (DEFAULT: GEO FILE NAME)
#             --parts N       ALSO WRITE N RANKS  (DEFAULT: NONE)
#             --clscale S     MULTIPLY ALL MESH SIZES BY S (2 = COARSER)
#             --threads T     GMSH MESHING THREADS (DEFAULT: OMP_NUM_THREADS OR 1)
#             --gmsh PATH     GMSH EXECUTABLE     (DEFAULT: THE gmsh NEXT TO
#                             --python, ELSE gmsh ON PATH. UBUNTU'S apt GMSH
#                             HAS NO NETGEN, SO Algorithm3D = 4 FAILS THERE.)
#             --python PATH   PYTHON WITH THE gmsh MODULE (DEFAULT: python3)
#             --overwrite     REPLACE AN EXISTING OUT_DIR
#
# EXAMPLE   : bash tools/cbs_case.sh cases/box/box.geo cases/box/box.par \
#                 cases/box/box.matprop ~/PhD/runs/box --parts 4
# ================================

set -euo pipefail

# ================================
# STEP 1: READ ARGUMENTS
# ================================
if [ "$#" -lt 4 ]; then
    sed -n '2,32p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

GEO="$(readlink -f "$1")"
PAR="$(readlink -f "$2")"
MATPROP="$(readlink -f "$3")"
OUT="$4"
shift 4

NAME="$(basename "$GEO")"
NAME="${NAME%.geo}"
NAME="${NAME%.msh}"
PARTS=0
CLSCALE=1
THREADS="${OMP_NUM_THREADS:-1}"
GMSH=""
PYTHON=python3
OVERWRITE=0

while [ "$#" -gt 0 ]; do
    case "$1" in
        --name)      NAME="$2";    shift 2 ;;
        --parts)     PARTS="$2";   shift 2 ;;
        --clscale)   CLSCALE="$2"; shift 2 ;;
        --threads)   THREADS="$2"; shift 2 ;;
        --gmsh)      GMSH="$2";    shift 2 ;;
        --python)    PYTHON="$2";  shift 2 ;;
        --overwrite) OVERWRITE=1;  shift ;;
        *) echo "ERROR: UNKNOWN OPTION $1" >&2; exit 1 ;;
    esac
done

TOOLS="$(cd "$(dirname "$0")" && pwd)"

# PREFER THE PIP GMSH THAT BELONGS TO THE CHOSEN PYTHON (FULL BUILD, SAME
# VERSION AS THE gmsh MODULE USED BY THE CONVERTER).
if [ -z "$GMSH" ]; then
    PYBIN="$(dirname "$(command -v "$PYTHON")")"
    if [ -x "$PYBIN/gmsh" ]; then
        GMSH="$PYBIN/gmsh"
    else
        GMSH=gmsh
    fi
fi

for f in "$GEO" "$PAR" "$MATPROP"; do
    [ -s "$f" ] || { echo "ERROR: MISSING OR EMPTY $f" >&2; exit 1; }
done
command -v "$GMSH" > /dev/null || { echo "ERROR: GMSH NOT FOUND: $GMSH" >&2; exit 1; }
"$PYTHON" -c "import gmsh" 2> /dev/null || {
    echo "ERROR: $PYTHON HAS NO gmsh MODULE (pip install gmsh, OR USE --python ~/venvs/phd/bin/python)" >&2
    exit 1
}

if [ -e "$OUT" ]; then
    [ "$OVERWRITE" -eq 1 ] || { echo "ERROR: $OUT EXISTS; USE --overwrite" >&2; exit 1; }
    rm -rf "$OUT"
fi
mkdir -p "$OUT/mesh" "$OUT/serial"
OUT="$(readlink -f "$OUT")"

# ================================
# STEP 2: GMSH .geo -> .msh
# ================================
# -3 MESHES THE VOLUME AND WRITES -o. (DO NOT ADD -save HERE: WITH A .geo IT
# WRITES THE FILE BEFORE MESHING, GIVING AN EMPTY MESH.)
if [[ "$GEO" == *.msh ]]; then
    echo "== 1/4 GMSH  : READY MESH $GEO (NOT RE-MESHED)"
    cp "$GEO" "$OUT/mesh/$NAME.msh"
    : > "$OUT/mesh/gmsh.log"
else
echo "== 1/4 GMSH  : $GEO"
"$GMSH" "$GEO" -3 -format msh2 -clscale "$CLSCALE" -nt "$THREADS" \
    -o "$OUT/mesh/$NAME.msh" > "$OUT/mesh/gmsh.log" 2>&1 || {
    tail -20 "$OUT/mesh/gmsh.log"; echo "ERROR: GMSH FAILED; SEE $OUT/mesh/gmsh.log" >&2; exit 1; }
fi
# GMSH CAN REPORT ERRORS (E.G. "No elements in volume") AND STILL EXIT 0.
if grep -q "^Error" "$OUT/mesh/gmsh.log"; then
    grep "^Error" "$OUT/mesh/gmsh.log" | head -5
    echo "ERROR: GMSH REPORTED ERRORS; SEE $OUT/mesh/gmsh.log" >&2
    exit 1
fi
echo "   gmsh $("$GMSH" --version 2>&1 | tail -1) ($GMSH)"
grep -E "^Info *: [0-9]+ nodes [0-9]+ elements" "$OUT/mesh/gmsh.log" | tail -1 || true

# ================================
# STEP 3: .msh -> SERIAL CBS3D CASE
# ================================
echo "== 2/4 CBS   : $OUT/serial/$NAME"
"$PYTHON" "$TOOLS/gmsh_to_cbs.py" "$OUT/mesh/$NAME.msh" \
    --case-name "$NAME" --out-dir "$OUT/serial" \
    --par "$PAR" --matprop "$MATPROP" | grep -v "^Info" || true
[ -s "$OUT/serial/$NAME.plt" ] || { echo "ERROR: CONVERSION FAILED" >&2; exit 1; }

# ================================
# STEP 4: VALIDATE
# ================================
echo "== 3/4 CHECK : mesh_summary.py"
"$PYTHON" "$TOOLS/mesh_summary.py" "$OUT/serial/$NAME" --json "$OUT/serial/mesh_summary.json"

# ================================
# STEP 5: OPTIONAL MPI PARTITIONS
# ================================
if [ "$PARTS" -gt 1 ]; then
    PDIR="$OUT/partitions_$PARTS"
    echo "== 4/4 MPI   : $PARTS RANKS -> $PDIR"
    "$PYTHON" "$TOOLS/export_rank_local_cbs.py" "$OUT/serial/$NAME" \
        --parts "$PARTS" --out-dir "$PDIR" --gmsh "$GMSH" | tail -3
    "$PYTHON" "$TOOLS/dd4d_check_partition.py" "$PDIR" "$NAME" "$PARTS" | tail -1
else
    echo "== 4/4 MPI   : SKIPPED (NO --parts)"
fi

# ================================
# STEP 6: HOW TO RUN
# ================================
echo
echo "CASE READY: $OUT"
echo "  SERIAL / OPENMP / PETSc (solver_opt in .par):"
echo "    cd $OUT/serial && OMP_NUM_THREADS=8 <build>/cbs3dpp_si $NAME"
if [ "$PARTS" -gt 1 ]; then
    echo "  MPI ($PARTS RANKS):"
    echo "    mkdir -p $OUT/run_p$PARTS && cd $OUT/run_p$PARTS"
    echo "    OMP_NUM_THREADS=1 mpirun -np $PARTS <mpi-build>/cbs3dpp_si $NAME $OUT/partitions_$PARTS"
fi
