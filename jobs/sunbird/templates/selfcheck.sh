#!/bin/bash
# ================================
# FILE      : selfcheck.sh
# PURPOSE   : CHECK THE CBS3D BUILD FOR MISTAKES BEFORE A PRODUCTION RUN:
#               1. THE TWO PRODUCTION EXECUTABLES EXIST AND MATCH THE SOURCE
#                  (BUILT AFTER THE NEWEST SOURCE/CMAKE FILE)
#               2. ALL UNIT / VERIFICATION TESTS (ctest, SEPARATE TEST BUILD)
#               3. SUBMIT selfcheck_slab.slurm: SERIAL + 4-RANK RUN OF A CASE
#                  WITH AN EXACT SOLUTION (NUCLEAR HEATING + WALL FLUX)
#
# USAGE     : bash $T/selfcheck.sh [JOBS]        (LOGIN NODE, ~3 MINUTES
#                                                 + A 4-CORE JOB OF ~10 MINUTES)
#
# RESULT    : "TESTS PASS" HERE, THEN "SELF-CHECK PASS" AT THE END OF
#             cbs3d_selfcheck_<JOBID>.out
# ================================

set -euo pipefail
source "$HOME/cbs3d_env.sh"
JOBS="${1:-8}"
T="$CBS3D_CODE/jobs/sunbird/templates"

# ================================
# STEP 1: EXECUTABLES NEWER THAN THE BUILD INPUTS
# ================================
# A COMMIT THAT TOUCHES ONLY DOCS OR JOB SCRIPTS DOES NOT CHANGE THE
# EXECUTABLES. COMPARE AGAINST THE NEWEST SOURCE/CMAKE FILE, NOT THE LAST
# COMMIT TIME: COMMIT f2ffbcc CHANGED ONLY selfcheck_slab.slurm AND MADE THE
# OLD CHECK REJECT A BUILD THAT WAS ALREADY CURRENT.
BUILD_INPUTS=(
    "$CBS3D_CODE/CMakeLists.txt"
    "$CBS3D_CODE/CMakePresets.json"
    "$CBS3D_CODE/apps"
    "$CBS3D_CODE/include"
    "$CBS3D_CODE/src"
)
LAST=$(find "${BUILD_INPUTS[@]}" -type f -printf '%T@\n' 2>/dev/null | sort -n | tail -1)
if [[ -z "$LAST" ]]; then
    echo "ERROR: no build inputs found under $CBS3D_CODE" >&2
    exit 1
fi
LAST=${LAST%.*}
for exe in "$CBS3D_SERIAL_EXE" "$CBS3D_MPI_EXE"; do
    [[ -x "$exe" ]] || { echo "ERROR: $exe MISSING; RUN build_cbs3d.sh" >&2; exit 1; }
    if [[ $(stat -c %Y "$exe") -lt "$LAST" ]]; then
        echo "ERROR: $exe IS OLDER THAN THE NEWEST BUILD INPUT ($(date -d "@$LAST" -Is)); RUN build_cbs3d.sh" >&2
        exit 1
    fi
done
echo "EXE    : both executables newer than the newest build input ($(date -d "@$LAST" -Is))"

# ================================
# STEP 2: UNIT AND VERIFICATION TESTS
# ================================
B="$CBS3D_CODE/build/selfcheck-tests"
export CC=mpicc CXX=mpicxx
if cmake -S "$CBS3D_CODE" -B "$B" -DCMAKE_BUILD_TYPE=Release \
        -DCBS3D_ENABLE_MPI=OFF -DCBS3D_ENABLE_PETSC=OFF -DCBS3D_ENABLE_OPENMP=OFF \
        -DCBS3D_BUILD_TESTS=ON -DCBS3D_BUILD_SOURCE_VALIDATION=ON > "$B.log" 2>&1 &&
   cmake --build "$B" --parallel "$JOBS" >> "$B.log" 2>&1 &&
   ( cd "$B" && ctest --output-on-failure ) >> "$B.log" 2>&1; then
    echo "TESTS  : PASS ($(grep -oE '[0-9]+ tests passed|tests passed, [0-9]+ tests failed out of [0-9]+' "$B.log" | tail -1))"
else
    tail -40 "$B.log"
    echo "TESTS  : FAIL (SEE $B.log)" >&2
    exit 1
fi

# ================================
# STEP 3: EXACT-SOLUTION RUNS ON A COMPUTE NODE
# ================================
JOB=$(sbatch --parsable "$T/selfcheck_slab.slurm")
echo "SLAB   : JOB $JOB (SERIAL + 4 RANKS AGAINST THE EXACT SOLUTION)"
echo "WAIT   : squeue -j $JOB ; then  tail -3 cbs3d_selfcheck_$JOB.out   (EXPECT: SELF-CHECK PASS)"
