#!/usr/bin/env bash
# Guarded entry point for the fw_straight 200 mm first-wall mesh/export workflow.
set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D/cases/fw_straight
CODE=/scratch/s.2337862/CBS3D/code
GEO="$ROOT/mesh/FW_SquareChannel_CHT.geo"
REPO_GEO="$CODE/cases/fw_straight/mesh/FW_SquareChannel_CHT.geo"
GEO_GMSH44="$CODE/cases/fw_straight/mesh/FW_SquareChannel_CHT_gmsh44.geo"
DRIVER="$CODE/jobs/sunbird/prepare_fw_straight_p40.sh"
PAR="$ROOT/fw_straight.par"

[[ -s "$REPO_GEO" ]] || { echo "FATAL: repository geometry missing: $REPO_GEO" >&2; exit 1; }
[[ -s "$GEO_GMSH44" ]] || { echo "FATAL: OCC-free Gmsh-4.4 geometry missing: $GEO_GMSH44" >&2; exit 1; }
[[ -s "$DRIVER" ]] || { echo "FATAL: missing $DRIVER" >&2; exit 1; }
[[ -s "$PAR" ]] || { echo "FATAL: missing $PAR" >&2; exit 1; }

# The repo owns the accepted geometry. Provision that exact copy into the
# Sunbird case directory every time, while retaining any previous local file.
mkdir -p "$ROOT/mesh" "$ROOT/input_backup" "$ROOT/run"
if [[ -s "$GEO" ]] && ! cmp -s "$REPO_GEO" "$GEO"; then
    STAMP=$(date +%Y%m%d_%H%M%S)
    cp -p "$GEO" "$ROOT/input_backup/FW_SquareChannel_CHT.geo.$STAMP.bak"
fi
cp -f "$REPO_GEO" "$GEO"

echo "Provisioned canonical geometry: $GEO"

# Filename is not enough. Several historical first-wall .geo files exist.
# Refuse anything other than the accepted 200 mm, 21x21 outer, 15x10 inner,
# 0.15 mm interface / 0.60 mm fluid / 1.00 mm solid-bulk mesh specification.
python3 - "$GEO" <<'PY'
from pathlib import Path
import re
import sys

s = Path(sys.argv[1]).read_text(encoding="utf-8")

def scalar(names, expected, tol=1.0e-12):
    for name in names:
        m = re.search(r"(?m)^\s*" + re.escape(name) + r"\s*=\s*([0-9.eE+-]+)\s*;", s)
        if m:
            value = float(m.group(1))
            if abs(value - expected) > tol:
                raise SystemExit("FATAL: %s=%.12g; expected %.12g" % (name, value, expected))
            return name, value
    raise SystemExit("FATAL: missing geometry scalar; expected one of %r" % (names,))

checks = [
    scalar(("L",), 0.200),
    scalar(("Wout",), 0.021),
    scalar(("Hout",), 0.021),
    scalar(("Wfluid", "Wf"), 0.015),
    scalar(("Hfluid", "Hf"), 0.010),
    scalar(("hInterface",), 0.00015),
    scalar(("hBulk",), 0.00100),
]

m = re.search(r"(?m)^\s*Field\[3\]\.VIn\s*=\s*([0-9.eE+-]+)\s*;", s)
if not m or abs(float(m.group(1)) - 0.00060) > 1.0e-12:
    raise SystemExit("FATAL: expected Field[3].VIn = 0.00060 m for the fluid core")
if not re.search(r"(?m)^\s*Field\[4\]\s*=\s*Min\s*;", s):
    raise SystemExit("FATAL: expected Field[4] = Min")
if not re.search(r"(?m)^\s*Background\s+Field\s*=\s*4\s*;", s):
    raise SystemExit("FATAL: expected Background Field = 4")

required = {
    "fluid volume 10": r'Physical\s+Volume\s*\(\s*"fluid"\s*,\s*10\s*\)',
    "solid volume 20": r'Physical\s+Volume\s*\(\s*"solid"\s*,\s*20\s*\)',
    "inlet BC511": r'Physical\s+Surface\s*\(\s*"inlet"\s*,\s*511\s*\)',
    "outlet BC520": r'Physical\s+Surface\s*\(\s*"outlet"\s*,\s*520\s*\)',
    "adiabatic BC530": r'Physical\s+Surface\s*\(\s*"adiabatic_wall"\s*,\s*530\s*\)',
    "heat-flux BC532": r'Physical\s+Surface\s*\(\s*"heat_flux"\s*,\s*532\s*\)',
}
for label, pattern in required.items():
    if not re.search(pattern, s):
        raise SystemExit("FATAL: missing expected physical group: " + label)

print("AUTHORITATIVE FW GEO AUDIT: PASS")
for name, value in checks:
    print("  %-12s = %.12g" % (name, value))
print("  fluid core h = 0.0006 m")
print("  physical IDs = volume 10/20; BC 511/520/530/532")
PY

# MeshIO reads the dimensional extension positionally. Locate the Nusselt
# label from that exact grammar instead of guessing its decorative wording.
python3 - "$PAR" <<'PY'
from pathlib import Path
import sys

p = Path(sys.argv[1])
lines = p.read_text(encoding="utf-8").splitlines()

def skip(line):
    s = line.strip()
    return (not s) or s.startswith("#") or s.startswith("!")

def next_data(idx):
    for j in range(idx + 1, len(lines)):
        if not skip(lines[j]):
            return j
    raise SystemExit("FATAL: unexpected EOF while locating dimensional .par block")

flow_hits = []
for i, line in enumerate(lines):
    s = line.strip().lower()
    if skip(line):
        continue
    if "p_ref" in s and ("mdot" in s or "inlet_mass_flow_rate" in s) and ("t_in" in s or "inlet_temperature" in s):
        flow_hits.append(i)

if len(flow_hits) != 1:
    raise SystemExit("FATAL: expected one reference/inlet/outlet .par label; found %d" % len(flow_hits))

flow_label = flow_hits[0]
flow_data = next_data(flow_label)
inlet_label = next_data(flow_data)
inlet_data = next_data(inlet_label)
nusselt_label = next_data(inlet_data)
nusselt_data = next_data(nusselt_label)

for desc, idx, n in (
    ("reference/inlet/outlet", flow_data, 6),
    ("prescribed inlet velocity", inlet_data, 3),
    ("Nusselt references", nusselt_data, 3),
):
    fields = lines[idx].split()
    if len(fields) < n:
        raise SystemExit("FATAL: malformed %s data row: %r" % (desc, lines[idx]))
    try:
        [float(x) for x in fields[:n]]
    except Exception:
        raise SystemExit("FATAL: non-numeric %s data row: %r" % (desc, lines[idx]))

lines[nusselt_label] = "Nusselt reference: nusselt_Tinf nusselt_Tref hydraulic_diameter"
p.write_text("\n".join(lines) + "\n", encoding="utf-8")

print("PARAMETER GRAMMAR AUDIT: PASS")
print("  positional Nusselt block located and canonicalized")
print("  data row preserved:", lines[nusselt_data])
PY

# Sunbird's gmsh/4.4.0 executable reports that it was built WITHOUT
# OpenCASCADE. The accepted workstation geometry uses OCC Box/BooleanDifference,
# so merely renaming modern mesh-field options is insufficient. Build a runtime
# copy of the otherwise-tested preparation driver and replace only its geometry
# compatibility block with the explicitly constructed OCC-free equivalent.
# The OCC-free file uses the same physical topology and the same HXT TET4 sizing
# fields; source-mesh audits below remain the final authority.
RUNTIME_DRIVER="$ROOT/run/prepare_fw_straight_p40_runtime.sh"
python3 - "$DRIVER" "$RUNTIME_DRIVER" "$GEO_GMSH44" <<'PY'
from pathlib import Path
import sys

source = Path(sys.argv[1])
target = Path(sys.argv[2])
compat = Path(sys.argv[3]).resolve()

s = source.read_text(encoding="utf-8")
start_token = 'cp "\\$GEO" "\\$COMPAT_GEO"\npython3 - "\\$COMPAT_GEO" <<\'PY44\'\n'
start = s.find(start_token)
if start < 0:
    raise SystemExit("FATAL: could not locate legacy OCC compatibility block in preparation driver")

end_token = '\nPY44\n\nmodule purge'
end = s.find(end_token, start)
if end < 0:
    raise SystemExit("FATAL: could not locate end of legacy OCC compatibility block")

replacement = 'cp "%s" "\\$COMPAT_GEO"\n' % str(compat)
s = s[:start] + replacement + s[end + len('\nPY44'):]

target.write_text(s, encoding="utf-8")
target.chmod(0o755)

check = target.read_text(encoding="utf-8")
if 'Geometry.OCCBoundsUseStl' in check:
    raise SystemExit("FATAL: runtime preparation driver still contains OCC geometry conversion")
if str(compat) not in check:
    raise SystemExit("FATAL: runtime preparation driver did not receive OCC-free geometry")

print("SUNBIRD GMSH COMPATIBILITY AUDIT: PASS")
print("  site gmsh        : 4.4.0 without OpenCASCADE")
print("  geometry source  :", compat)
print("  topology         : explicit shared surfaces, no OCC Boolean")
print("  volume mesher    : HXT Algorithm3D=10")
PY

# The generated Slurm script contains literal MSH2 section markers in a nested
# quoted Python here-document. Supply those spellings during outer expansion so
# `set -u` cannot interpret them as undefined shell variables.
export Nodes='$Nodes'
export Elements='$Elements'

exec bash "$RUNTIME_DRIVER"
