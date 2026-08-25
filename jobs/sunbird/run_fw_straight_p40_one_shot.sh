#!/usr/bin/env bash
# Guarded entry point for the fw_straight 200 mm first-wall mesh/export workflow.
set -euo pipefail

ROOT=/scratch/s.2337862/CBS3D/cases/fw_straight
CODE=/scratch/s.2337862/CBS3D/code
GEO="$ROOT/mesh/FW_SquareChannel_CHT.geo"
REPO_GEO="$CODE/cases/fw_straight/mesh/FW_SquareChannel_CHT.geo"
DRIVER="$CODE/jobs/sunbird/prepare_fw_straight_p40.sh"

[[ -s "$REPO_GEO" ]] || { echo "FATAL: repository geometry missing: $REPO_GEO" >&2; exit 1; }
[[ -s "$DRIVER" ]] || { echo "FATAL: missing $DRIVER" >&2; exit 1; }

# The previous launcher incorrectly assumed the case-local .geo already existed.
# The repo now owns the accepted geometry. Provision that exact copy into the
# Sunbird case directory every time, while retaining any previous local file.
mkdir -p "$ROOT/mesh" "$ROOT/input_backup"
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

exec bash "$DRIVER"
