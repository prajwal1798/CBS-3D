#!/usr/bin/env python3
"""G8a: compare legacy-VTU and native-checkpoint restart on the SAME MPI mesh.

Independent owned-node global-ID comparison of pressure, velocity, temperature
at step 50000 (loaded state) and step 50001 (one coupled production CBS step).
No source numerics, restart files or mesh inputs are modified.
This tests native-vs-legacy two-rank restart parity, NOT rank-count invariance.
"""

import argparse
import math
from pathlib import Path
import xml.etree.ElementTree as ET

from audit_temperature_undershoot_vtu import read_piece, files_for_step, named, floats


def fields(directory, step):
    rank_files = files_for_step(directory, step)
    if len(rank_files) != 2:
        raise ValueError(f"{directory}: expected 2 MPI VTU pieces at step {step}")
    owned = {}
    all_gids = set()
    for file in rank_files:
        piece = read_piece(file)
        root = ET.parse(file).getroot()
        pdata = root.find("./UnstructuredGrid/Piece/PointData")
        pressure = floats(named(pdata, ("pressure",), str(file)), str(file))
        if len(pressure) != len(piece["gids"]):
            raise ValueError(f"{file}: invalid pressure length")
        for i, gid in enumerate(piece["gids"]):
            all_gids.add(gid)
            if piece["owned"][i] != 1:
                continue
            if gid in owned:
                raise ValueError(f"{directory}: duplicated owned node {gid}")
            owned[gid] = {
                "temperature": (piece["temperature"][i],),
                "velocity": piece["velocity"][i],
                "pressure": (pressure[i],),
                "coords": piece["points"][i],
            }
    if set(owned) != all_gids:
        raise ValueError("VTU owned nodes do not cover all observed IDs")
    return owned


def compare(reference, candidate, step, tolerance):
    a = fields(reference, step)
    b = fields(candidate, step)
    if a.keys() != b.keys():
        missing_a = b.keys()-a.keys()
        missing_b = a.keys()-b.keys()
        raise ValueError(f"global ID mismatch step={step}: extra={len(missing_a)}, missing={len(missing_b)}")

    errors = {}
    for name in ("temperature", "velocity", "pressure", "coords"):
        best = (0.0, -1)
        for gid in a:
            for x, y in zip(a[gid][name], b[gid][name]):
                diff = abs(x-y)
                if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(diff)):
                    raise ValueError(f"nonfinite {name} at gid={gid}")
                if diff > best[0]:
                    best = (diff, gid)
        errors[name] = best
    print(f"step {step}: {len(a)} unique owned nodes")
    for name, (diff, gid) in errors.items():
        print(f"  max_abs_{name:11s}={diff:.12e}, worst_gid={gid}")
        if diff > tolerance[name]:
            raise ValueError(
                f"G8a PARITY FAIL at step {step}, {name}: "
                f"{diff:.12e} > {tolerance[name]:.12e}"
            )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--legacy-vtu-dir", type=Path, required=True)
    ap.add_argument("--native-vtu-dir", type=Path, required=True)
    ap.add_argument("--before", type=int, default=50000)
    ap.add_argument("--after", type=int, default=50001)
    ap.add_argument("--temperature-tol", type=float, default=1e-9)
    ap.add_argument("--velocity-tol", type=float, default=1e-10)
    ap.add_argument("--pressure-tol", type=float, default=1e-7)
    args = ap.parse_args()
    tolerance = {
        "temperature": args.temperature_tol,
        "velocity": args.velocity_tol,
        "pressure": args.pressure_tol,
        "coords": 1e-12,
    }
    print("=== G8a native vs legacy restart, 2-rank coupled solver ===")
    print(f"Legacy VTU root: {args.legacy_vtu_dir}")
    print(f"Native VTU root: {args.native_vtu_dir}")
    for step in (args.before, args.after):
        compare(args.legacy_vtu_dir, args.native_vtu_dir, step, tolerance)
    print("G8A NATIVE / LEGACY RESTART PARITY PASS")
    print("NOTE: This is same-mesh, same-2-rank restart parity; serial/partition-count parity is not tested.")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, IndexError) as exc:
        raise SystemExit(f"G8a PARITY AUDIT ERROR: {exc}")
