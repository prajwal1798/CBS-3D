#!/bin/bash
# ================================
# FILE      : blanket.sh
# PURPOSE   : EVERYTHING FOR THE LAMINAR BLANKET-SLICE STUDY ON SUNBIRD2,
#             ONE SUB-COMMAND PER STEP (LOGIN NODE; HEAVY WORK GOES TO SLURM):
#
#   update [FILE ...]   FILE = A GIT BUNDLE (.bundle: ITS BRANCH IS CHECKED OUT
#                       AS claude/blanket) OR PATCH FILES (git am); WITHOUT
#                       ARGUMENTS git pull --ff-only. THEN REBUILD AND RUN
#                       selfcheck.sh
#   setup MSH [--force] CONVERT THE GMSH MESH TO A CBS CASE (ONE JOB) AND
#                       BUILD THE PARTITION LADDER $LADDER (ONE JOB EACH,
#                       STARTED WHEN THE CONVERSION HAS FINISHED)
#   smoke [STAGE]       300 ITERATIONS ON THE SMALLEST SET (30 MIN LIMIT):
#                       START-UP, HEAT INPUT, TIME STEP, NO CRASH
#   run STAGE [OPTS]    STAGE = L0 (FLOW ONLY), L1 (+ WALL FLUX), L2 (+ NUCLEAR
#                       HEATING), L1T (L1 IN PHYSICAL TIME, FOR TIME AVERAGES).
#                       PICKS THE LARGEST READY SET THAT FITS THE FREE CORES
#                       (submit_auto.sh). OPTS: --time T --max N --sets N --dry-run
#   energy RUN_DIR [L1|L2]  IMPLICIT STEADY ENERGY (WATER + STEEL, ONE PETSc
#                       SOLVE, SteadyEnergySolver) ON THE LAST FLOW FIELD OF
#                       RUN_DIR (ITS mean_velocity WHEN IT HAS ONE), SAME RANKS
#   check [RUN_DIR]     MASS + ENERGY BALANCE AND RESIDUAL PLOT OF A RUN
#                       (DEFAULT: THE NEWEST RUN)
#   status              YOUR JOBS AND THE LATEST PROGRESS LINE OF EACH RUN
#
# SETTINGS  : ENVIRONMENT VARIABLES (DEFAULTS IN BRACKETS)
#   BLANKET_DIR   CASE FOLDER          [$CBS3D_SCRATCH/cases/blanket_slice]
#   BLANKET_CASE  CASE NAME            [blanket]
#   LADDER        PARTITION SETS       ["40 80 128 256"]
#   CP_SCALE      DIVIDE THE SOLID c_p BY THIS FACTOR [1]. THE STEADY RESULT
#                 DOES NOT DEPEND ON THE SOLID c_p; A FACTOR > 1 SHORTENS THE
#                 SOLID'S THERMAL TIME CONSTANT (TRANSIENT IS THEN NOT PHYSICAL)
#   CBS3D_LTS_MONOTONE     [1] LOCAL TIME STEPS MAY ONLY DECREASE ON A REFRESH
#                 (OTHERWISE THE 3.2 M-TET RUN REFRESHED EVERY ~2.5 ITERATIONS
#                 AND CYCLED INSTEAD OF CONVERGING)
#   CBS3D_AVERAGE_START_TIME  [OFF] TRANSIENT STAGE L1T: TIME AVERAGES FROM
#                 THIS PHYSICAL TIME [s] (TimeAverage.hpp); check THEN ALSO
#                 BALANCES THE MEAN FIELDS
#   CBS3D_TEMPERATURE_AFC  [1] ALGEBRAIC FLUX CORRECTION OF THE TEMPERATURE
#                 (NO UNDER/OVERSHOOT AT HIGH CELL PECLET NUMBER); 0 = OFF
#
# EXAMPLE   : bash $T/blanket.sh update ~/cbs3d-blanket.bundle
#             bash $T/blanket.sh setup $CBS3D_SCRATCH/meshes/slice_h5.msh
#             bash $T/blanket.sh smoke L1
#             bash $T/blanket.sh run L1 --time 48:00:00
#             bash $T/blanket.sh check
# ================================

set -euo pipefail
usage() { sed -n '2,48p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
[[ "$#" -ge 1 ]] || usage

source "$HOME/cbs3d_env.sh"
T="$CBS3D_CODE/jobs/sunbird/templates"
SRC="$CBS3D_CODE/cases/blanket_slice"
DIR="${BLANKET_DIR:-$CBS3D_SCRATCH/cases/blanket_slice}"
CASE="${BLANKET_CASE:-blanket}"
LADDER="${LADDER:-40 80 128 256}"
CP_SCALE="${CP_SCALE:-1}"
# WATER AT Pr ~ 6 GIVES CELL PECLET NUMBERS OF SEVERAL HUNDRED: WITHOUT THE
# LIMITER THE TEMPERATURE UNDERSHOOTS THE INLET VALUE BY TENS OF KELVIN.
export CBS3D_TEMPERATURE_AFC="${CBS3D_TEMPERATURE_AFC:-1}"
export CBS3D_LTS_MONOTONE="${CBS3D_LTS_MONOTONE:-1}"

CMD="$1"
shift

# ================================
# HELPERS
# ================================
stage_par() {
    local f="$SRC/blanket_$1.par"
    [[ -s "$f" ]] || { echo "ERROR: UNKNOWN STAGE '$1' (NO $f); USE L0, L1, L1T OR L2" >&2; exit 1; }
    echo "$f"
}

# MATPROP FOR THE RUN: THE CASE COPY, OR ONE WITH EVERY SOLID c_p / CP_SCALE.
run_matprop() {
    if [[ "$CP_SCALE" == 1 ]]; then
        echo "$SRC/blanket.matprop"
        return
    fi
    local out="$DIR/blanket_cpscale${CP_SCALE}.matprop"
    awk -v s="$CP_SCALE" 'NR == 1 || /^#/ { print; next }
        $3 == "solid" { $5 = $5 / s } { print }' "$SRC/blanket.matprop" > "$out"
    echo "$out"
}

newest_run() {
    ls -1dt "$DIR"/run_"${CASE}"_* 2>/dev/null | head -1
}

progress() {
    # LAST "it=" PROGRESS LINE OF A solver.out, WITHOUT COLOURS
    sed 's/\x1b\[[0-9;]*m//g' "$1" 2>/dev/null | tr '\r' '\n' |
        grep -oE "(it=|Iteration )[0-9]+/[0-9]+.*" | tail -1 | cut -c1-170
}

case "$CMD" in

# ================================
# UPDATE THE CODE, BUILD, SELF-CHECK
# ================================
update)
    FILES=()
    for f in "$@"; do FILES+=("$(readlink -f "$f")"); done
    set -- "${FILES[@]}"
    cd "$CBS3D_CODE"
    if [[ -n "$(git status --porcelain --untracked-files=no)" ]]; then
        echo "ERROR: $CBS3D_CODE HAS UNCOMMITTED CHANGES; COMMIT OR STASH THEM FIRST" >&2
        git status --short --untracked-files=no
        exit 1
    fi
    if [[ "$#" -eq 1 && "$1" == *.bundle ]]; then
        B="$1"
        git bundle verify "$B" > /dev/null || { echo "ERROR: $B IS NOT A VALID BUNDLE FOR THIS REPOSITORY" >&2; exit 1; }
        HEAD_REF="$(git bundle list-heads "$B" | awk 'NR == 1 { print $2 }')"
        git fetch "$B" "+$HEAD_REF:refs/heads/claude/blanket"
        git checkout claude/blanket
        echo "BRANCH : claude/blanket FROM $B ($HEAD_REF); YOUR PREVIOUS BRANCH IS UNCHANGED"
    elif [[ "$#" -gt 0 ]]; then
        for p in "$@"; do
            [[ -s "$p" ]] || { echo "ERROR: NO PATCH FILE $p" >&2; exit 1; }
        done
        git config user.email > /dev/null || {
            echo "ERROR: git NEEDS YOUR IDENTITY TO APPLY PATCHES, ONCE:" >&2
            echo "  git config --global user.name  \"Your Name\"" >&2
            echo "  git config --global user.email \"you@example.com\"" >&2
            exit 1; }
        git am --3way "$@" || {
            echo "ERROR: PATCH DID NOT APPLY. UNDO WITH: git -C $CBS3D_CODE am --abort" >&2
            exit 1; }
    else
        git pull --ff-only
    fi
    git log --oneline -3
    bash "$T/build_cbs3d.sh"
    bash "$T/selfcheck.sh"
    ;;

# ================================
# MESH -> CBS CASE -> PARTITION LADDER
# ================================
setup)
    MSH="$(readlink -f "${1:?USAGE: blanket.sh setup MSH [--force]}")"
    FORCE="${2:-}"
    [[ -s "$MSH" ]] || { echo "ERROR: NO MESH $MSH" >&2; exit 1; }
    if compgen -G "$DIR/partitions_*" > /dev/null || compgen -G "$DIR/run_*" > /dev/null; then
        [[ "$FORCE" == --force ]] || {
            echo "ERROR: $DIR ALREADY HAS PARTITIONS OR RUNS. setup WOULD DELETE THEM." >&2
            echo "       USE ANOTHER BLANKET_DIR, OR ADD --force." >&2
            exit 1; }
    fi
    mkdir -p "$DIR"
    PARENT="$(dirname "$DIR")"
    # THE CONVERSION EMPTIES AND RECREATES $DIR, SO ITS LOG GOES ONE LEVEL UP.
    J=$(cd "$PARENT" && sbatch --parsable \
        --output="$PARENT/cbs3d_case_%j.out" --error="$PARENT/cbs3d_case_%j.err" \
        "$T/build_case.slurm" "$MSH" "$(stage_par L1)" \
        "$SRC/blanket.matprop" "$DIR" --name "$CASE" --overwrite)
    echo "CONVERT: JOB $J  $MSH -> $DIR/serial/$CASE  (LOG $PARENT/cbs3d_case_$J.out)"
    # THE PARTITION JOBS START ONLY AFTER THE CONVERSION HAS SUCCEEDED.
    bash "$T/prepare_partitions.sh" "$DIR" "$CASE" $LADDER --after "$J"
    echo "NEXT   : squeue -u \$USER ; WHEN ALL ARE DONE: bash $T/blanket.sh smoke L1"
    ;;

# ================================
# SHORT START-UP TEST
# ================================
smoke)
    STAGE="${1:-L1}"
    PAR="$DIR/smoke_${STAGE}.par"
    # 300 ITERATIONS, NO EARLY STOP, ONE FINAL VTU
    awk '/^#/ { print; next }
         want == 1 { $1 = 300; want = 0 }
         want == 2 { $NF = 1000000; $6 = 300; want = 0 }
         /ntime transient_on/ { want = 1 }
         /steady_min_iterations/ { want = 2 }
         { print }' "$(stage_par "$STAGE")" > "$PAR"
    SMALL="$(echo $LADDER | tr ' ' '\n' | sort -n | head -1)"
    bash "$T/submit_auto.sh" "$DIR" "$CASE" --sets "$SMALL" --time 00:30:00 \
        --par "$PAR" --matprop "$(run_matprop)"
    echo "CHECK  : bash $T/blanket.sh status   (THEN LOOK FOR HEAT_INPUT_TOTAL_W, dt AND NO ERROR)"
    ;;

# ================================
# PRODUCTION RUN
# ================================
run)
    STAGE="${1:?USAGE: blanket.sh run L0|L1|L2 [--time T] [--max N] [--dry-run]}"
    shift
    PAR="$(stage_par "$STAGE")"
    MAT="$(run_matprop)"
    echo "STAGE  : $STAGE  ($PAR)"
    echo "SOLID  : c_p divided by $CP_SCALE  ($MAT)"
    echo "AFC    : CBS3D_TEMPERATURE_AFC=$CBS3D_TEMPERATURE_AFC   LTS: CBS3D_LTS_MONOTONE=$CBS3D_LTS_MONOTONE"
    echo "AVERAGE: start time ${CBS3D_AVERAGE_START_TIME:-off} s, start iteration ${CBS3D_AVERAGE_START_ITERATION:-off}"
    if [[ "${CBS3D_RESTART:-0}" != 0 ]]; then
        echo "RESTART: ${CBS3D_RESTART_FORMAT:-native} FROM ${CBS3D_RESTART_ROOT:-?} (ITERATION ${CBS3D_RESTART_ITERATION:-from manifest}); AS INITIAL CONDITION: ${CBS3D_RESTART_AS_INITIAL_CONDITION:-0}"
        echo "         THE RUN MUST USE THE SAME RANK COUNT (PASS --sets N)"
    fi
    bash "$T/submit_auto.sh" "$DIR" "$CASE" --time 24:00:00 "$@" --par "$PAR" --matprop "$MAT"
    ;;

# ================================
# IMPLICIT STEADY ENERGY ON THE FLOW OF AN EARLIER RUN
# ================================
energy)
    RUN="$(readlink -f "${1:?USAGE: blanket.sh energy RUN_DIR [L1|L2] [OPTS]}")"
    STAGE="${2:-L1}"
    shift; [[ "$#" -gt 0 ]] && shift
    LAST="$(ls -1 "$RUN"/output/"$CASE"/"$CASE"_step_*.pvtu 2>/dev/null | sort | tail -1)"
    [[ -n "$LAST" ]] || { echo "ERROR: NO .pvtu IN $RUN/output/$CASE" >&2; exit 1; }
    ITER="$(basename "$LAST" .pvtu | sed 's/.*_step_0*//')"
    ARR=velocity
    grep -q 'Name="mean_velocity"' "$LAST" && ARR=mean_velocity
    NR=$(ls -d "$RUN"/partitions/rank_* | wc -l)
    export CBS3D_STEADY_ENERGY=1 CBS3D_RESTART=1 CBS3D_RESTART_FORMAT=legacy_vtu
    export CBS3D_RESTART_ROOT="$RUN/output/$CASE" CBS3D_RESTART_ITERATION="$ITER"
    export CBS3D_RESTART_VELOCITY_ARRAY="$ARR"
    unset CBS3D_AVERAGE_START_TIME CBS3D_AVERAGE_START_ITERATION
    echo "FLOW   : $LAST ($ARR, ITERATION $ITER, $NR RANKS)"
    echo "STAGE  : $STAGE  ($(stage_par "$STAGE")); AFC $CBS3D_TEMPERATURE_AFC"
    bash "$T/submit_auto.sh" "$DIR" "$CASE" --sets "$NR" --time 01:00:00 "$@" \
        --par "$(stage_par "$STAGE")" --matprop "$(run_matprop)"
    echo "THEN   : bash $T/blanket.sh status ; bash $T/blanket.sh check"
    ;;

# ================================
# ENERGY BALANCE + RESIDUALS
# ================================
check)
    RUN="${1:-$(newest_run)}"
    [[ -n "$RUN" && -d "$RUN" ]] || { echo "ERROR: NO RUN FOLDER IN $DIR" >&2; exit 1; }
    RUN="$(readlink -f "$RUN")"
    echo "RUN    : $RUN"
    progress "$RUN/solver.out" || true
    sed 's/\x1b\[[0-9;]*m//g' "$RUN/solver.out" | grep -E "HEAT INPUT|nuclear heating|BC 532|HEAT_INPUT_TOTAL_W|stop reason" || true
    if ls "$RUN"/output/"$CASE"/*.pvtu > /dev/null 2>&1; then
        python3 "$CBS3D_CODE/tools/energy_balance.py" "$RUN/output/$CASE" \
            --log "$RUN/solver.out" --plane y=0 --rho 997 --cp 4182 || true
        # TRANSIENT RUN WITH TIME AVERAGES: ALSO THE BALANCE OF THE MEAN FIELDS
        LAST="$(ls -1t "$RUN"/output/"$CASE"/*.pvtu | head -1)"
        if grep -q 'Name="mean_velocity"' "$LAST"; then
            echo "---- TIME-AVERAGED FIELDS ----"
            sed 's/\x1b\[[0-9;]*m//g' "$RUN/solver.out" | grep "\[AVG\]" | tail -1 || true
            python3 "$CBS3D_CODE/tools/energy_balance.py" "$LAST" --mean \
                --log "$RUN/solver.out" --plane y=0 --rho 997 --cp 4182 || true
        fi
    else
        echo "NOTE   : NO .pvtu YET (FIRST ONE AT THE vtu_output_every_iterations STEP)"
    fi
    if ls "$RUN"/output/"$CASE"/*residuals.csv > /dev/null 2>&1; then
        python3 "$CBS3D_CODE/tools/plot_residuals_pub.py" "$RUN/output/$CASE" \
            --out "$RUN/residuals" && echo "PLOT   : $RUN/residuals.png"
    fi
    ;;

# ================================
# WHAT IS RUNNING
# ================================
status)
    squeue -u "$USER" -o "%.10i %.16j %.8T %.10M %.6D %.6C %R" || true
    for r in $(ls -1dt "$DIR"/run_"${CASE}"_* 2>/dev/null | head -5); do
        echo "$(basename "$r"): $(progress "$r/solver.out")"
        grep -h "ERROR\|FATAL" "$r/solver.out" "$r/solver.err" 2>/dev/null | head -3 || true
    done
    ;;

*)
    usage
    ;;
esac
