#!/bin/bash
# ==============================================================================
# One-shot Sunbird launcher for the ONLY active turbulence verification case:
# wall-resolved standard-SA NASA/TMR zero-pressure-gradient flat plate.
#
# It rebuilds the current integration/sa-mpi solver, audits the retained P0040
# partition as a wall-resolved y+~1 benchmark, and submits the production run.
# If the retained partition is actually a coarse wall-model mesh, this script
# refuses to submit rather than produce an impressively precise wrong answer.
# ==============================================================================

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D
CODE=$ROOT/code
CASE=$ROOT/cases/flatplate
PART=$CASE/partition/coarse/P0040/cbs
EXE=$CODE/build/mpi-release/cbs3d_parallel
JOB=$CODE/jobs/sunbird/flatplate_sa_tmr_p40.slurm

GCC_ROOT=/lustrehome/apps/compilers/gnu/12.1.0
MPI_ROOT=/apps/libraries/openmpi/4.1.6/el7/AVX512/gnu-12.1
PETSC_DIR=/home/s.2337862/software/petsc-gcc12-install

mkdir -p "$CASE/jobs/logs"

printf '===== ACTIVE SOURCE =====\n'
cd "$CODE"
printf 'branch : %s\n' "$(git rev-parse --abbrev-ref HEAD)"
printf 'HEAD   : %s\n' "$(git rev-parse HEAD)"
[[ "$(git rev-parse --abbrev-ref HEAD)" == "integration/sa-mpi" ]] || {
    echo "FATAL: expected integration/sa-mpi" >&2
    exit 1
}

if [[ -n "$(git status --porcelain)" ]]; then
    echo "FATAL: worktree is dirty after synchronisation; refusing benchmark build." >&2
    git status --short
    exit 1
fi

printf '\n===== BUILD MPI/PETSC RELEASE =====\n'
module purge
module load cmake/3.31.10
export PATH="$MPI_ROOT/bin:$GCC_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$PETSC_DIR/lib:$MPI_ROOT/lib:$GCC_ROOT/lib64:${LD_LIBRARY_PATH:-}"
export PETSC_DIR
unset PETSC_ARCH
export OMPI_CC="$GCC_ROOT/bin/gcc"
export OMPI_CXX="$GCC_ROOT/bin/g++"

cmake --preset mpi-release
cmake --build --preset build-mpi-release --parallel 4
[[ -x "$EXE" ]] || { echo "FATAL: executable missing after build: $EXE" >&2; exit 1; }
ls -lh "$EXE"

printf '\n===== P40 STRUCTURE + PHYSICS + WALL-RESOLUTION AUDIT =====\n'
python3 - "$PART" <<'PY'
from pathlib import Path
import math
import statistics
import sys

root = Path(sys.argv[1])
if not root.is_dir():
    raise SystemExit("FATAL: partition directory missing: {}".format(root))

ranks = sorted(p for p in root.glob("rank_[0-9][0-9][0-9][0-9]") if p.is_dir())
if len(ranks) != 40:
    raise SystemExit("FATAL: expected 40 rank directories, found {}".format(len(ranks)))

for ext in ("plt", "bco", "mpi", "par", "matprop"):
    count = sum(1 for r in ranks if len(list(r.glob("*." + ext))) == 1)
    if count != 40:
        raise SystemExit("FATAL: expected exactly one .{} per rank; valid ranks={}".format(ext, count))

rep = ranks[0]
par_path = next(rep.glob("*.par"))
mat_path = next(rep.glob("*.matprop"))

lines = par_path.read_text(encoding="utf-8").splitlines()

def next_data(label_pred, description):
    hits = [i for i, s in enumerate(lines) if label_pred(s.strip().lower())]
    if len(hits) != 1:
        raise SystemExit("FATAL: expected one {} label in {}; found {}".format(description, par_path, len(hits)))
    for j in range(hits[0] + 1, len(lines)):
        t = lines[j].strip()
        if t and not t.startswith(("#", "!")):
            return t.split()
    raise SystemExit("FATAL: missing data after {}".format(description))

free = next_data(lambda s: "freestream" in s and "u" in s, "freestream")
reline = next_data(lambda s: s.startswith("re ") or "re pr ra ri" in s, "Re/Pr")
dim = next_data(lambda s: "dimensional" in s and "mass-flow" in s, "dimensional controls")
sa = next_data(lambda s: "spalart-allmaras controls" in s, "SA controls")
sain = next_data(lambda s: "spalart-allmaras inlet" in s, "SA inlet controls")
art = next_data(lambda s: "artificial diffusion" in s, "artificial diffusion")

if len(free) < 5 or abs(float(free[0]) - 1.0) > 1e-12 or any(abs(float(free[k])) > 1e-12 for k in (1, 2)):
    raise SystemExit("FATAL: TMR benchmark requires Uinf=(1,0,0); got {}".format(free[:3]))
if len(reline) < 1 or abs(float(reline[0]) - 5.0e6) > 1.0:
    raise SystemExit("FATAL: TMR benchmark requires Re_x(x=1)=5e6; got {}".format(reline[0] if reline else "missing"))
if dim[:3] != ["1", "1", "0"]:
    raise SystemExit("FATAL: expected dimensional/material/mass-flow flags [1 1 0]; got {}".format(dim[:3]))
if sa[:3] != ["1", "0", "0"]:
    raise SystemExit("FATAL: expected standard SA [1 0 0]; got {}".format(sa[:3]))
if not sain or abs(float(sain[0]) - 3.0) > 1e-12:
    raise SystemExit("FATAL: expected nu_tilde_inf/nu=3; got {}".format(sain[0] if sain else "missing"))
if not art or int(float(art[0])) != 0:
    raise SystemExit("FATAL: TMR benchmark requires art_diff=0")

mat = [s.strip() for s in mat_path.read_text(encoding="utf-8").splitlines()
       if s.strip() and not s.lstrip().startswith(("#", "!"))]
if len(mat) < 2 or int(mat[0].split()[0]) != 1:
    raise SystemExit("FATAL: expected one all-fluid material in {}".format(mat_path))
f = mat[1].split()
if len(f) < 8 or int(f[0]) != 0 or f[2].lower() != "fluid":
    raise SystemExit("FATAL: malformed fluid .matprop row: {}".format(mat[1]))
rho = float(f[3]); mu = float(f[6]); nu = mu / rho
if abs(rho - 1.0) > 1e-12 or abs(nu - 2.0e-7) > 2.0e-15:
    raise SystemExit("FATAL: TMR scaling requires rho=1 and nu=2e-7; got rho={} nu={}".format(rho, nu))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def norm(a):
    return math.sqrt(dot(a,a))

def percentile(values, p):
    x = sorted(values)
    if not x:
        return float("nan")
    q = (len(x)-1)*p/100.0
    i = int(math.floor(q)); j = min(i+1, len(x)-1); t = q-i
    return x[i]*(1-t) + x[j]*t

wall_yplus = []
wall_h = []
bc_seen = set()
wall_faces = 0

for rank in ranks:
    bco_path = next(rank.glob("*.bco"))
    blines = [s.strip() for s in bco_path.read_text(encoding="utf-8").splitlines()
              if s.strip() and not s.lstrip().startswith(("#", "!"))]
    nflag = int(blines[0].split()[0])
    mapping = {}
    for row in blines[1:1+nflag]:
        p = row.split(); mapping[int(p[0])] = int(p[1])
    bc_seen.update(mapping.values())

    plt_path = next(rank.glob("*.plt"))
    with plt_path.open("r", encoding="utf-8") as fh:
        header = fh.readline().split()
        if len(header) < 3:
            raise SystemExit("FATAL: bad .plt header in {}".format(plt_path))
        nelem, npoin, nboun = map(int, header[:3])
        tets = [None]*(nelem+1)
        for _ in range(nelem):
            p = fh.readline().split(); ie = int(p[0]); tets[ie] = tuple(map(int,p[1:5]))
        xyz = [None]*(npoin+1)
        for _ in range(npoin):
            p = fh.readline().split(); ip = int(p[0]); xyz[ip] = tuple(map(float,p[1:4]))
        for _ in range(nboun):
            p = fh.readline().split()
            face = tuple(map(int,p[:3])); parent = int(p[3]); raw = int(p[4])
            if mapping.get(raw) != 530:
                continue
            pts = [xyz[n] for n in face]
            xc = sum(v[0] for v in pts)/3.0
            if xc <= 0.0:
                continue
            tet = tets[parent]
            tp = [xyz[n] for n in tet]
            detj = abs(dot(sub(tp[1],tp[0]), cross(sub(tp[2],tp[0]), sub(tp[3],tp[0]))))
            area = 0.5*norm(cross(sub(pts[1],pts[0]), sub(pts[2],pts[0])))
            if not (detj > 0.0 and area > 0.0):
                raise SystemExit("FATAL: degenerate wall-adjacent tetrahedron in {}".format(plt_path))
            h = detj/(2.0*area)
            wall_faces += 1
            if 0.8 <= xc <= 1.2:
                rex = 5.0e6*xc
                cf = 0.0592*rex**(-0.2)
                utau = math.sqrt(0.5*cf)
                wall_h.append(h)
                wall_yplus.append(h*utau/nu)

for bc in (506, 511, 520, 530):
    if bc not in bc_seen:
        raise SystemExit("FATAL: global P40 BC inventory missing {}".format(bc))
if wall_faces <= 0 or len(wall_yplus) < 4:
    raise SystemExit("FATAL: could not establish wall-resolution statistics near x=1")

median_h = statistics.median(wall_h)
median_yp = statistics.median(wall_yplus)
p95_yp = percentile(wall_yplus,95)
max_yp = max(wall_yplus)

print("rank directories               = 40")
print("U_inf                          = 1")
print("Re_x(x=1)                      = 5e6")
print("nu                             = {:.12e}".format(nu))
print("nu_tilde_inf/nu                = 3")
print("wall model                     = OFF")
print("physical BC inventory          = {}".format(sorted(bc_seen)))
print("BC530 physical wall faces      = {}".format(wall_faces))
print("x=0.8..1.2 wall samples        = {}".format(len(wall_yplus)))
print("median first-TET altitude      = {:.12e}".format(median_h))
print("design y+ median near x=1      = {:.6f}".format(median_yp))
print("design y+ P95 near x=1         = {:.6f}".format(p95_yp))
print("design y+ max near x=1         = {:.6f}".format(max_yp))

# Wall-resolved SA benchmark gate.  The estimate uses only the standard design
# correlation; the accepted y+ will be recomputed from converged wall shear.
if median_yp > 2.0 or p95_yp > 3.0:
    raise SystemExit(
        "FATAL: retained P0040 is not a y+~1 wall-resolved TMR mesh "
        "(median={:.3f}, P95={:.3f}). Do not benchmark SA on it.".format(median_yp,p95_yp)
    )

print("P40 WALL-RESOLVED TMR PREFLIGHT: PASS")
PY

printf '\n===== SUBMIT BENCHMARK =====\n'
JOBID=$(sbatch --parsable "$JOB")
printf 'JOB ID  : %s\n' "$JOBID"
printf 'OUTPUT  : %s/jobs/logs/sa_tmr_%s.out\n' "$CASE" "$JOBID"
printf 'ERROR   : %s/jobs/logs/sa_tmr_%s.err\n' "$CASE" "$JOBID"
printf 'RUN ROOT: %s/runs/sa_tmr_wall_resolved/job_%s\n' "$CASE" "$JOBID"
squeue -j "$JOBID" || true
