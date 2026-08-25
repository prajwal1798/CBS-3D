#!/usr/bin/env bash
# ==============================================================================
# CBS3D++_SI | Sunbird one-shot first-wall mesh/input/P40 exporter
#
# Case:
#   /scratch/s.2337862/CBS3D/cases/fw_straight
#
# Authoritative geometry:
#   mesh/FW_SquareChannel_CHT.geo
#
# This driver:
#   1. preserves the user's existing .par numerical controls and material values;
#   2. patches only first-wall physics/CHT/SA fields that are known and required;
#   3. makes the solid volumetric source unambiguous (.par only, .matprop Qvol=0);
#   4. computes Re, Pr and Ubulk from the retained helium material properties;
#   5. submits one legal 40-core Sunbird Slurm job;
#   6. generates Gmsh-4.4-compatible source MSH2;
#   7. partitions with METIS while preserving source MSH2 element numbering;
#   8. uses the established pure-Python MSH2 exporter to write 40 CBS rank trees;
#   9. hard-audits tags, material IDs, rank files and total owned TET4 count.
#
# IMPORTANT:
#   - It does NOT invent helium rho/cp/k/mu. Those are retained from the existing
#     authoritative fw_straight.matprop supplied for this engineering state.
#   - It does NOT run the solver. It prepares the immutable P0040 CBS inputs.
# ==============================================================================

set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D/cases/fw_straight
CODE=/scratch/s.2337862/CBS3D/code
MESH_DIR="$ROOT/mesh"
GEO="$MESH_DIR/FW_SquareChannel_CHT.geo"
PAR="$ROOT/fw_straight.par"
MATPROP="$ROOT/fw_straight.matprop"
P40="$ROOT/partition/fine/P0040"
RUN_DIR="$ROOT/run"
CASE=fw_straight
NPART=40

for f in "$GEO" "$PAR" "$MATPROP"; do
    [[ -s "$f" ]] || { echo "FATAL: missing/empty required input: $f" >&2; exit 1; }
done

mkdir -p "$ROOT/input_backup" "$RUN_DIR"
STAMP=$(date +%Y%m%d_%H%M%S)
cp -p "$PAR" "$ROOT/input_backup/fw_straight.par.$STAMP.bak"
cp -p "$MATPROP" "$ROOT/input_backup/fw_straight.matprop.$STAMP.bak"

# ----------------------------------------------------------------------------
# Patch and validate the authoritative .par/.matprop on the login node.
# Keep the existing parser order and all unrelated numerical controls intact.
# ----------------------------------------------------------------------------
python3 - "$PAR" "$MATPROP" <<'PY'
from __future__ import print_function
from pathlib import Path
import math
import re
import sys

par_path = Path(sys.argv[1])
mat_path = Path(sys.argv[2])

MDOT = 0.05
TIN = 573.15
PREF = 8.0e6
Q_SOLID = 1.093e7
Q_FLUX = 5.0e5
AREA = 0.015 * 0.010
DH = 2.0 * 0.015 * 0.010 / (0.015 + 0.010)
EUROFER_K = 28.3
PR_T = 0.90
SA_INLET_RATIO = 3.0

raw = mat_path.read_text(encoding="utf-8").splitlines()
data_idx = [i for i, line in enumerate(raw)
            if line.strip() and not line.lstrip().startswith(("#", "!"))]
if len(data_idx) < 3:
    raise SystemExit("FATAL: fw_straight.matprop must contain nmat + two material rows")

try:
    nmat = int(raw[data_idx[0]].split()[0])
except Exception:
    raise SystemExit("FATAL: invalid material-count line in fw_straight.matprop")
if nmat != 2:
    raise SystemExit("FATAL: expected exactly 2 materials (helium=0, EUROFER=1), found %d" % nmat)

rows = {}
row_indices = {}
for idx in data_idx[1:1+nmat]:
    fields = raw[idx].split()
    if len(fields) < 8:
        raise SystemExit(
            "FATAL: .matprop row must begin: id name phase rho cp k mu Qvol; got %d fields in %r"
            % (len(fields), raw[idx]))
    try:
        mid = int(fields[0])
        rho, cp, k, mu, qvol = map(float, fields[3:8])
    except Exception:
        raise SystemExit("FATAL: invalid numeric material row: %r" % raw[idx])
    if mid in rows:
        raise SystemExit("FATAL: duplicate material ID %d" % mid)
    if rho <= 0.0 or cp <= 0.0 or k <= 0.0 or mu < 0.0:
        raise SystemExit("FATAL: material %d has invalid rho/cp/k/mu" % mid)
    rows[mid] = [fields[0], fields[1], fields[2], rho, cp, k, mu, qvol]
    row_indices[mid] = idx

if set(rows) != {0, 1}:
    raise SystemExit("FATAL: material IDs must be exactly {0,1}; found %r" % sorted(rows))
if rows[0][2].lower() != "fluid":
    raise SystemExit("FATAL: material 0 must be phase 'fluid'")
if rows[1][2].lower() != "solid":
    raise SystemExit("FATAL: material 1 must be phase 'solid'")

rho_h, cp_h, k_h, mu_h = rows[0][3:7]
if mu_h <= 0.0:
    raise SystemExit("FATAL: helium viscosity must be positive for Reynolds/SA calculation")

rho_s, cp_s, _k_s, mu_s = rows[1][3:7]
if abs(mu_s) > 1.0e-30:
    raise SystemExit("FATAL: solid viscosity must be zero; found %.12g" % mu_s)

def fmt_mat(row, rho, cp, k, mu, qvol):
    return "%s %s %s %.12g %.12g %.12g %.12g %.12g" % (
        row[0], row[1], row[2], rho, cp, k, mu, qvol)

raw[row_indices[0]] = fmt_mat(rows[0], rho_h, cp_h, k_h, mu_h, 0.0)
raw[row_indices[1]] = fmt_mat(rows[1], rho_s, cp_s, EUROFER_K, 0.0, 0.0)
mat_path.write_text("\n".join(raw) + "\n", encoding="utf-8")

UBULK = MDOT / (rho_h * AREA)
RE = MDOT * DH / (AREA * mu_h)
PR = mu_h * cp_h / k_h

if not all(math.isfinite(v) and v > 0.0 for v in (UBULK, RE, PR)):
    raise SystemExit("FATAL: non-finite derived helium Ubulk/Re/Pr")

lines = par_path.read_text(encoding="utf-8").splitlines()

def is_comment_or_blank(s):
    t = s.strip()
    return (not t) or t.startswith("#") or t.startswith("!")

def data_line_after(label_idx):
    for j in range(label_idx + 1, len(lines)):
        if not is_comment_or_blank(lines[j]):
            return j
    raise SystemExit("FATAL: no data line after .par label: %s" % lines[label_idx])

def find_label(predicate, description):
    hits = [i for i, s in enumerate(lines) if not is_comment_or_blank(s) and predicate(s.strip())]
    if len(hits) != 1:
        raise SystemExit("FATAL: expected one .par label for %s; found %d" % (description, len(hits)))
    return hits[0]

def replace_vector(label_pred, description, values, expected_min=None):
    i = find_label(label_pred, description)
    j = data_line_after(i)
    old = lines[j].split()
    if expected_min is not None and len(old) < expected_min:
        raise SystemExit("FATAL: malformed data row for %s" % description)
    lines[j] = "  ".join(str(v) for v in values)
    return old

def replace_first_token(label_pred, description, value):
    i = find_label(label_pred, description)
    j = data_line_after(i)
    fields = lines[j].split()
    if not fields:
        raise SystemExit("FATAL: empty data row for %s" % description)
    fields[0] = str(value)
    lines[j] = "  ".join(fields)

replace_first_token(
    lambda s: "TEMP_CALC" in s.upper() or "TEMP CALC" in s.upper(),
    "TEMP_CALC", 1)

replace_vector(
    lambda s: all(tok in s.upper() for tok in ("UX", "UY", "UZ")) and ("P" in s.upper()) and ("T" in s.upper()),
    "freestream Ux Uy Uz P T",
    ["%.12g" % UBULK, "0.0", "0.0", "0.0", "%.12g" % TIN],
    expected_min=5)

replace_vector(
    lambda s: re.search(r"\bRE\b", s, re.I) and re.search(r"\bPR\b", s, re.I)
              and re.search(r"\bRA\b", s, re.I) and re.search(r"\bRI\b", s, re.I),
    "Re Pr Ra Ri",
    ["%.12g" % RE, "%.12g" % PR, "0.0", "0.0"],
    expected_min=4)

replace_vector(
    lambda s: "source_solid" in s and "heat_flux_bc" in s,
    "alpha_sf k_ratio source_solid heat_flux_bc",
    ["1.0", "1.0", "%.12g" % Q_SOLID, "%.12g" % Q_FLUX],
    expected_min=4)

replace_vector(
    lambda s: "dimensional_mode" in s and "material_properties_enabled" in s and "mass_flow_inlet_enabled" in s,
    "dimensional/material/mass-flow controls",
    ["1", "1", "1"],
    expected_min=3)

i = find_label(
    lambda s: "p_ref" in s and ("inlet_mass_flow_rate" in s or "mdot" in s)
              and ("inlet_temperature" in s or "T_in" in s),
    "reference/inlet/outlet controls")
j = data_line_after(i)
flow = lines[j].split()
if len(flow) < 6:
    raise SystemExit("FATAL: malformed p_ref/model_depth/mdot/rho/T/outlet row")
model_depth = float(flow[1])
if not math.isfinite(model_depth) or model_depth <= 0.0:
    raise SystemExit("FATAL: model_depth must remain positive")
lines[j] = "  ".join([
    "%.12g" % PREF,
    "%.12g" % model_depth,
    "%.12g" % MDOT,
    "%.12g" % rho_h,
    "%.12g" % TIN,
    "0.0",
])

try:
    replace_vector(
        lambda s: "inlet_u" in s and "inlet_v" in s and "inlet_w" in s,
        "fallback inlet velocity",
        ["%.12g" % UBULK, "0.0", "0.0"],
        expected_min=3)
except SystemExit as exc:
    if "found 0" not in str(exc):
        raise
    replace_vector(
        lambda s: "prescribed inlet velocity" in s.lower(),
        "fallback inlet velocity",
        ["%.12g" % UBULK, "0.0", "0.0"],
        expected_min=3)

i = find_label(
    lambda s: "Nusselt" in s and ("diameter" in s.lower() or "reference" in s.lower()),
    "Nusselt reference")
j = data_line_after(i)
nu = lines[j].split()
if len(nu) < 3:
    raise SystemExit("FATAL: malformed Nusselt reference row")
nu[2] = "%.12g" % DH
lines[j] = "  ".join(nu)

replace_vector(
    lambda s: "turbulence_on" in s and "turbulence_model" in s and "turbulent_thermal_diffusivity_on" in s,
    "Spalart-Allmaras controls",
    ["1", "0", "1"],
    expected_min=3)

replace_vector(
    lambda s: "sa_inlet_ratio" in s and "sa_prandtl_t" in s,
    "Spalart-Allmaras inlet/thermal controls",
    ["%.12g" % SA_INLET_RATIO, "%.12g" % PR_T],
    expected_min=2)

replace_vector(
    lambda s: "sa_min_wall_distance" in s and "sa_min_stilde" in s and "sa_nu_tilde_floor" in s,
    "Spalart-Allmaras numerical floors",
    ["1.0e-14", "1.0e-14", "0.0"],
    expected_min=3)

replace_vector(
    lambda s: "sa_use_stilde_limiter" in s and "sa_implicit_destruction" in s,
    "Spalart-Allmaras switches",
    ["1", "1", "1.0e6"],
    expected_min=2)

par_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

print("FIRST-WALL INPUT PATCH: PASS")
print("  geometry fluid area     = %.12g m^2" % AREA)
print("  hydraulic diameter      = %.12g m" % DH)
print("  mdot                    = %.12g kg/s" % MDOT)
print("  T_in                    = %.12g K" % TIN)
print("  p_ref                   = %.12g Pa" % PREF)
print("  helium rho              = %.12g kg/m^3" % rho_h)
print("  helium cp               = %.12g J/(kg K)" % cp_h)
print("  helium k                = %.12g W/(m K)" % k_h)
print("  helium mu               = %.12g Pa s" % mu_h)
print("  U_bulk                  = %.12g m/s" % UBULK)
print("  Re(Dh)                  = %.12g" % RE)
print("  Pr                      = %.12g" % PR)
print("  EUROFER rho             = %.12g kg/m^3" % rho_s)
print("  EUROFER cp              = %.12g J/(kg K)" % cp_s)
print("  EUROFER k               = %.12g W/(m K)" % EUROFER_K)
print("  source_solid(.par)      = %.12g W/m^3" % Q_SOLID)
print("  heat_flux_bc(.par)      = %.12g W/m^2" % Q_FLUX)
print("  Qvol helium/.matprop    = 0")
print("  Qvol EUROFER/.matprop   = 0")
print("  SA / Pr_t               = standard / %.12g" % PR_T)
PY

python3 - "$PAR" "$MATPROP" <<'PY'
from pathlib import Path
import sys

par = Path(sys.argv[1]).read_text().splitlines()
mat = [x.strip() for x in Path(sys.argv[2]).read_text().splitlines()
       if x.strip() and not x.lstrip().startswith(("#", "!"))]

def next_values(token):
    for i, line in enumerate(par):
        if token in line:
            for j in range(i + 1, len(par)):
                t = par[j].strip()
                if t and not t.startswith(("#", "!")):
                    return t.split()
    raise SystemExit("FATAL: missing label containing %s" % token)

src = next_values("source_solid")
sa = next_values("turbulence_on")
dim = next_values("mass_flow_inlet_enabled")
flow = next_values("p_ref")

if len(mat) != 3:
    raise SystemExit("FATAL: expected 3 data lines in two-material .matprop")
q0 = float(mat[1].split()[7])
q1 = float(mat[2].split()[7])
if float(src[2]) != 1.093e7 or float(src[3]) != 5.0e5:
    raise SystemExit("FATAL: incorrect first-wall .par source/flux after patch")
if q0 != 0.0 or q1 != 0.0:
    raise SystemExit("FATAL: .matprop Qvol must be zero for both materials")
if sa[:3] != ["1", "0", "1"]:
    raise SystemExit("FATAL: expected standard SA + turbulent thermal diffusivity")
if dim[:3] != ["1", "1", "1"]:
    raise SystemExit("FATAL: expected dimensional/material/mass-flow flags = 1 1 1")
if abs(float(flow[0]) - 8.0e6) > 1.0 or abs(float(flow[2]) - 0.05) > 1e-14 or abs(float(flow[4]) - 573.15) > 1e-10:
    raise SystemExit("FATAL: p_ref/mdot/Tin audit failed")
print("POST-PATCH INPUT AUDIT: PASS")
PY

EXPORTER=$(find /scratch/s.2337862/CBS3D -type f -name 'export_rank_local_cbs_msh2_py36.py' -print -quit 2>/dev/null || true)
[[ -n "$EXPORTER" && -s "$EXPORTER" ]] || {
    echo "FATAL: export_rank_local_cbs_msh2_py36.py not found under /scratch/s.2337862/CBS3D" >&2
    exit 1
}

echo "Exporter: $EXPORTER"

JOBFILE="$RUN_DIR/prepare_fw_straight_p40_${STAMP}.slurm"

cat > "$JOBFILE" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=fw_prepare_p40
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --output=$RUN_DIR/fw_prepare_p40_%j.out
#SBATCH --error=$RUN_DIR/fw_prepare_p40_%j.err

set -euo pipefail

ROOT="$ROOT"
CODE="$CODE"
GEO="$GEO"
PAR="$PAR"
MATPROP="$MATPROP"
P40="$P40"
CASE="$CASE"
NPART="$NPART"
EXPORTER="$EXPORTER"

ARCHIVE="\$ROOT/partition/fine/archive"
GMSH_DIR="\$P40/gmsh"
CBS_DIR="\$P40/cbs"
LOG_DIR="\$P40/logs"
SOURCE_MSH="\$GMSH_DIR/\${CASE}_source.msh"
PART_MSH="\$GMSH_DIR/\${CASE}_partitioned_40.msh"
COMPAT_GEO="\$GMSH_DIR/FW_SquareChannel_CHT_sunbird_gmsh44.geo"

printf '============================================================\n'
printf 'CBS3D FW_STRAIGHT P40 PREPARE/EXPORT\n'
printf 'job       = %s\n' "\$SLURM_JOB_ID"
printf 'host      = %s\n' "\$(hostname)"
printf 'date      = %s\n' "\$(date -Is)"
printf 'geo       = %s\n' "\$GEO"
printf 'exporter  = %s\n' "\$EXPORTER"
printf 'target    = %s\n' "\$CBS_DIR"
printf '============================================================\n'

if [[ -e "\$P40" ]]; then
    mkdir -p "\$ARCHIVE"
    NOW=\$(date +%Y%m%d_%H%M%S)
    mv "\$P40" "\$ARCHIVE/P0040_\${NOW}"
fi
mkdir -p "\$GMSH_DIR" "\$CBS_DIR" "\$LOG_DIR"

cp "\$GEO" "\$COMPAT_GEO"
python3 - "\$COMPAT_GEO" <<'PY44'
from pathlib import Path
import re
import sys

p = Path(sys.argv[1])
s = p.read_text(encoding="utf-8")
replacements = {
    "Field[1].SurfacesList": "Field[1].FacesList",
    "Field[2].SizeMin": "Field[2].LcMin",
    "Field[2].SizeMax": "Field[2].LcMax",
    "Mesh.MeshSizeFromPoints": "Mesh.CharacteristicLengthFromPoints",
    "Mesh.MeshSizeFromCurvature": "Mesh.CharacteristicLengthFromCurvature",
    "Mesh.MeshSizeExtendFromBoundary": "Mesh.CharacteristicLengthExtendFromBoundary",
    "Mesh.MeshSizeMin": "Mesh.CharacteristicLengthMin",
    "Mesh.MeshSizeMax": "Mesh.CharacteristicLengthMax",
}
for old, new in replacements.items():
    s = s.replace(old, new)
s = re.sub(r"^\s*Field\[1\]\.Sampling\s*=\s*[^;]+;\s*$", "", s, flags=re.M)
p.write_text(s, encoding="utf-8")
PY44

module purge
module load gmsh/4.4.0
export OMP_NUM_THREADS=40

printf '\n===== GMSH VERSION =====\n'
gmsh -version

printf '\n===== SOURCE MESH =====\n'
gmsh "\$COMPAT_GEO" -3 -format msh2 -o "\$SOURCE_MSH" 2>&1 | tee "\$LOG_DIR/gmsh_generate.log"
[[ -s "\$SOURCE_MSH" ]] || { echo 'FATAL: source MSH2 missing/empty' >&2; exit 1; }

python3 - "\$SOURCE_MSH" <<'PYMESH' | tee "\$LOG_DIR/source_mesh_audit.txt"
from __future__ import print_function
import sys

path = sys.argv[1]
nodes = None
elements = None
tets = 0
tris = 0
other_volume = 0
by_tet = {}
by_tri = {}
with open(path, "r") as f:
    it = iter(f)
    for line in it:
        key = line.strip()
        if key == "$Nodes":
            nodes = int(next(it).strip())
        elif key == "$Elements":
            elements = int(next(it).strip())
            for _ in range(elements):
                p = next(it).split()
                etype = int(p[1]); ntags = int(p[2])
                phys = int(p[3]) if ntags >= 1 else 0
                if etype == 4:
                    tets += 1; by_tet[phys] = by_tet.get(phys, 0) + 1
                elif etype == 2:
                    tris += 1; by_tri[phys] = by_tri.get(phys, 0) + 1
                elif etype in (5, 6, 7, 11, 12, 17):
                    other_volume += 1
            break

print("nodes                =", nodes)
print("elements             =", elements)
print("TET4                 =", tets)
print("TRI3                 =", tris)
print("tet physical counts  =", sorted(by_tet.items()))
print("tri physical counts  =", sorted(by_tri.items()))
print("other volume elems   =", other_volume)
if not nodes or not elements or tets <= 0:
    raise SystemExit("FATAL: no usable TET4 source mesh")
if other_volume != 0:
    raise SystemExit("FATAL: non-TET4 volume elements detected")
for tag in (10, 20):
    if by_tet.get(tag, 0) <= 0:
        raise SystemExit("FATAL: missing volume physical tag %d" % tag)
for tag in (511, 520, 530, 532):
    if by_tri.get(tag, 0) <= 0:
        raise SystemExit("FATAL: missing boundary physical tag %d" % tag)
print("SOURCE MESH AUDIT: PASS")
print("SOURCE_TETS=%d" % tets)
PYMESH

SOURCE_TETS=\$(awk -F= '/^TET4/{gsub(/[[:space:]]/,"",\$2); print \$2}' "\$LOG_DIR/source_mesh_audit.txt")
[[ "\$SOURCE_TETS" =~ ^[0-9]+$ && "\$SOURCE_TETS" -gt 0 ]] || { echo 'FATAL: could not recover source TET4 count' >&2; exit 1; }

printf '\n===== METIS P40 PARTITION =====\n'
gmsh "\$SOURCE_MSH" \
    -part "\$NPART" \
    -format msh2 \
    -save_all \
    -preserve_numbering_msh2 \
    -setnumber Mesh.PartitionOldStyleMsh2 1 \
    -setnumber Mesh.PartitionCreateGhostCells 0 \
    -setnumber Mesh.PartitionSplitMeshFiles 0 \
    -o "\$PART_MSH" \
    2>&1 | tee "\$LOG_DIR/gmsh_partition.log"

[[ -s "\$PART_MSH" ]] || { echo 'FATAL: partitioned MSH2 missing/empty' >&2; exit 1; }
PART_ELEMS=\$(python3 -c 'import sys; f=open(sys.argv[1]); it=iter(f); target=chr(36)+"Elements"; print(next(it).strip()) if any(line.strip()==target for line in it) else None' "\$PART_MSH")
[[ "\$PART_ELEMS" =~ ^[0-9]+$ && "\$PART_ELEMS" -gt 0 ]] || { echo 'FATAL: partitioned MSH2 has zero/invalid element count' >&2; exit 1; }
echo "partitioned MSH2 elements = \$PART_ELEMS"

printf '\n===== CBS RANK EXPORT =====\n'
python3 - "\$EXPORTER" "\$SOURCE_MSH" "\$PART_MSH" "\$CBS_DIR" "\$PAR" "\$MATPROP" "\$CASE" <<'PYEXP'
from __future__ import print_function
import subprocess
import sys

exporter, source, partitioned, outdir, par, matprop, case = sys.argv[1:]
py = sys.executable
help_text = subprocess.check_output([py, exporter, "--help"], stderr=subprocess.STDOUT, universal_newlines=True)
cmd = [py, exporter]

def choose(flags, value, label):
    for flag in flags:
        if flag in help_text:
            cmd.extend([flag, value]); return True
    raise SystemExit("FATAL: exporter option not recognized for %s; inspect --help" % label)

source_flag = next((f for f in ("--source-msh", "--source_msh", "--source-mesh", "--source") if f in help_text), None)
part_flag = next((f for f in ("--partitioned-msh", "--partitioned_msh", "--partitioned-mesh", "--partition-msh", "--partitioned") if f in help_text), None)
if bool(source_flag) != bool(part_flag):
    raise SystemExit("FATAL: exporter exposes only one of source/partitioned mesh options")
if source_flag:
    cmd += [source_flag, source, part_flag, partitioned]
else:
    cmd += [source, partitioned]

choose(("--parts", "--nparts", "--partitions"), "40", "partition count")
choose(("--case-name", "--case_name", "--case", "--name"), case, "case name")
choose(("--out-dir", "--out_dir", "--output-dir", "--output"), outdir, "output directory")
choose(("--par",), par, ".par")
choose(("--matprop",), matprop, ".matprop")

print("RUNNING:")
print(" ".join(cmd))
subprocess.check_call(cmd)
PYEXP

printf '\n===== HARD P40 CBS AUDIT =====\n'
NRANK=\$(find "\$CBS_DIR" -maxdepth 1 -type d -name 'rank_[0-9][0-9][0-9][0-9]' | wc -l)
[[ "\$NRANK" -eq 40 ]] || { echo "FATAL: expected 40 rank directories, found \$NRANK" >&2; exit 1; }

SUM_TETS=0
for i in \$(seq 0 39); do
    r=\$(printf '%04d' "\$i")
    d="\$CBS_DIR/rank_\$r"
    for ext in plt bco material mpi par matprop; do
        f="\$d/\${CASE}_rank_\${r}.\${ext}"
        [[ -s "\$f" ]] || { echo "FATAL: missing/empty \$f" >&2; exit 1; }
    done
    nt=\$(awk 'NR==1{print \$1; exit}' "\$d/\${CASE}_rank_\${r}.plt")
    [[ "\$nt" =~ ^[0-9]+$ ]] || { echo "FATAL: invalid .plt header on rank \$r" >&2; exit 1; }
    SUM_TETS=\$((SUM_TETS + nt))
done

[[ "\$SUM_TETS" -eq "\$SOURCE_TETS" ]] || {
    echo "FATAL: sum(rank owned TET4)=\$SUM_TETS != source TET4=\$SOURCE_TETS" >&2
    exit 1
}

PAR_HASHES=\$(find "\$CBS_DIR" -name '*.par' -exec sha256sum {} \; | awk '{print \$1}' | sort -u | wc -l)
MAT_HASHES=\$(find "\$CBS_DIR" -name '*.matprop' -exec sha256sum {} \; | awk '{print \$1}' | sort -u | wc -l)
[[ "\$PAR_HASHES" -eq 1 ]] || { echo 'FATAL: rank .par files differ' >&2; exit 1; }
[[ "\$MAT_HASHES" -eq 1 ]] || { echo 'FATAL: rank .matprop files differ' >&2; exit 1; }

[[ "\$(sha256sum "\$PAR" | awk '{print \$1}')" == "\$(sha256sum "\$CBS_DIR/rank_0000/\${CASE}_rank_0000.par" | awk '{print \$1}')" ]] || {
    echo 'FATAL: exported .par does not match authoritative root .par' >&2; exit 1; }
[[ "\$(sha256sum "\$MATPROP" | awk '{print \$1}')" == "\$(sha256sum "\$CBS_DIR/rank_0000/\${CASE}_rank_0000.matprop" | awk '{print \$1}')" ]] || {
    echo 'FATAL: exported .matprop does not match authoritative root .matprop' >&2; exit 1; }

MATS=\$(awk 'NF>=6 {print \$NF}' "\$CBS_DIR"/rank_*/*.material | sort -n -u | tr '\n' ' ')
echo "material IDs across P40 = \$MATS"
echo "\$MATS" | grep -Eq '(^| )0( |\$)' || { echo 'FATAL: fluid material ID 0 absent' >&2; exit 1; }
echo "\$MATS" | grep -Eq '(^| )1( |\$)' || { echo 'FATAL: solid material ID 1 absent' >&2; exit 1; }

BCS=\$(awk 'NF>=2 && \$1 ~ /^[0-9]+\$/ && \$2 ~ /^[0-9]+\$/ {print \$2}' "\$CBS_DIR"/rank_*/*.bco | sort -n -u | tr '\n' ' ')
echo "solver BC IDs across P40 = \$BCS"
for bc in 511 520 530 532; do
    echo "\$BCS" | grep -Eq "(^| )\${bc}( |\$)" || { echo "FATAL: BC \$bc absent from exported P40" >&2; exit 1; }
done

echo "rank directories          = \$NRANK"
echo "source TET4               = \$SOURCE_TETS"
echo "sum owned rank TET4       = \$SUM_TETS"
echo "CBS root                  = \$CBS_DIR"
echo
printf '============================================================\n'
printf 'FW_STRAIGHT P40 MESH + CBS EXPORT: PASS\n'
printf 'Permanent CBS inputs: %s\n' "\$CBS_DIR"
printf '============================================================\n'
EOF

JOBID=$(sbatch --parsable "$JOBFILE")

echo
echo "============================================================"
echo "SUBMITTED ONE-SHOT FW_STRAIGHT P40 PREPARATION"
echo "JOB ID : $JOBID"
echo "JOBFILE: $JOBFILE"
echo "OUTPUT : $RUN_DIR/fw_prepare_p40_${JOBID}.out"
echo "ERROR  : $RUN_DIR/fw_prepare_p40_${JOBID}.err"
echo "TARGET : $P40/cbs"
echo "============================================================"
squeue -j "$JOBID" || true
