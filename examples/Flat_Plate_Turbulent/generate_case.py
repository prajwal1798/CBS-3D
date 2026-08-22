#!/usr/bin/env python3
"""Generate a 3-D thin-span TET4 zero-pressure-gradient flat-plate case.

The benchmark follows the NASA TMR flat-plate scaling used for SA verification:
Re_x = 5e6 at x=1 and nu_tilde_inf/nu_inf = 3.  The mesh is a structured
Cartesian brick lattice split with a globally consistent six-tetrahedron
Freudenthal decomposition.  The first wall-normal interval is derived from a
requested target y+ using Cf = 0.0592 Re_x^(-1/5) at x_ref=1.

The script writes solver-native .plt, .bco, .par, .material and .matprop files.
No Gmsh conversion or external Python package is required.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

Point = Tuple[float, float, float]
Tet = Tuple[int, int, int, int]
Face = Tuple[int, int, int]


def geometric_ratio(first: float, total: float, intervals: int) -> float:
    if first <= 0.0 or total <= 0.0 or intervals < 1:
        raise ValueError("invalid geometric spacing inputs")
    if intervals * first >= total:
        return 1.0

    def series(r: float) -> float:
        if abs(r - 1.0) < 1.0e-14:
            return intervals * first
        return first * (r**intervals - 1.0) / (r - 1.0)

    lo, hi = 1.0, 2.0
    while series(hi) < total:
        hi *= 1.5
        if hi > 100.0:
            raise RuntimeError("failed to bracket wall-normal growth ratio")

    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if series(mid) < total:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def wall_normal_coordinates(first: float, height: float, intervals: int) -> Tuple[List[float], float]:
    ratio = geometric_ratio(first, height, intervals)
    coords = [0.0]
    dy = first
    for _ in range(intervals):
        coords.append(coords[-1] + dy)
        dy *= ratio
    coords[-1] = height
    return coords, ratio


def streamwise_coordinates(xmin: float, xmax: float, nup: int, nplate: int) -> List[float]:
    if not (xmin < 0.0 < xmax):
        raise ValueError("flat-plate domain must straddle x=0")
    if nup < 1 or nplate < 2:
        raise ValueError("insufficient streamwise intervals")

    upstream = [xmin + (0.0 - xmin) * i / nup for i in range(nup + 1)]
    exponent = 1.5
    plate = [xmax * (i / nplate) ** exponent for i in range(1, nplate + 1)]
    return upstream + plate


def node_id(i: int, j: int, k: int, nx: int, ny: int) -> int:
    return 1 + i + (nx + 1) * (j + (ny + 1) * k)


def make_mesh(xs: Sequence[float], ys: Sequence[float], zs: Sequence[float]) -> Tuple[List[Point], List[Tet]]:
    nx = len(xs) - 1
    ny = len(ys) - 1
    nz = len(zs) - 1

    points: List[Point] = []
    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                expected = node_id(i, j, k, nx, ny)
                if expected != len(points) + 1:
                    raise RuntimeError("node numbering inconsistency")
                points.append((x, y, z))

    tets: List[Tet] = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                v000 = node_id(i, j, k, nx, ny)
                v100 = node_id(i + 1, j, k, nx, ny)
                v010 = node_id(i, j + 1, k, nx, ny)
                v110 = node_id(i + 1, j + 1, k, nx, ny)
                v001 = node_id(i, j, k + 1, nx, ny)
                v101 = node_id(i + 1, j, k + 1, nx, ny)
                v011 = node_id(i, j + 1, k + 1, nx, ny)
                v111 = node_id(i + 1, j + 1, k + 1, nx, ny)

                tets.extend(
                    [
                        (v000, v100, v110, v111),
                        (v000, v110, v010, v111),
                        (v000, v010, v011, v111),
                        (v000, v011, v001, v111),
                        (v000, v001, v101, v111),
                        (v000, v101, v100, v111),
                    ]
                )

    return points, tets


def signed_det(points: Sequence[Point], tet: Tet) -> float:
    p = [points[n - 1] for n in tet]
    ax, ay, az = (p[1][q] - p[0][q] for q in range(3))
    bx, by, bz = (p[2][q] - p[0][q] for q in range(3))
    cx, cy, cz = (p[3][q] - p[0][q] for q in range(3))
    return (
        ax * (by * cz - bz * cy)
        - bx * (ay * cz - az * cy)
        + cx * (ay * bz - az * by)
    )


def boundary_faces(points: Sequence[Point], tets: Sequence[Tet]) -> List[Tuple[Face, int, int]]:
    inventory: Dict[Tuple[int, int, int], List[object]] = {}
    local_faces = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))

    for ie, tet in enumerate(tets, start=1):
        for local in local_faces:
            face = tuple(tet[q] for q in local)
            key = tuple(sorted(face))
            if key not in inventory:
                inventory[key] = [face, ie, 1]
            else:
                inventory[key][2] = int(inventory[key][2]) + 1

    xmin = min(p[0] for p in points)
    xmax = max(p[0] for p in points)
    ymin = min(p[1] for p in points)
    ymax = max(p[1] for p in points)
    zmin = min(p[2] for p in points)
    zmax = max(p[2] for p in points)
    scale = max(xmax - xmin, ymax - ymin, zmax - zmin, 1.0)
    tol = 2.0e-11 * scale

    result: List[Tuple[Face, int, int]] = []
    for face, parent, multiplicity in inventory.values():
        if multiplicity != 1:
            continue

        f = tuple(int(v) for v in face)
        xyz = [points[n - 1] for n in f]
        cx = sum(p[0] for p in xyz) / 3.0

        if all(abs(p[0] - xmin) <= tol for p in xyz):
            bc = 511
        elif all(abs(p[0] - xmax) <= tol for p in xyz):
            bc = 520
        elif all(abs(p[1] - ymax) <= tol for p in xyz):
            bc = 506
        elif all(abs(p[1] - ymin) <= tol for p in xyz):
            bc = 506 if cx < -tol else 530
        elif all(abs(p[2] - zmin) <= tol for p in xyz) or all(
            abs(p[2] - zmax) <= tol for p in xyz
        ):
            bc = 506
        else:
            raise RuntimeError(f"unclassified boundary face {f}")

        result.append((f, int(parent), bc))

    return result


def write_plt(path: Path, points: Sequence[Point], tets: Sequence[Tet], faces: Sequence[Tuple[Face, int, int]]) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write(f"{len(tets)} {len(points)} {len(faces)}\n")
        for ie, tet in enumerate(tets, start=1):
            out.write(f"{ie} {tet[0]} {tet[1]} {tet[2]} {tet[3]}\n")
        for ip, (x, y, z) in enumerate(points, start=1):
            out.write(f"{ip} {x:.17e} {y:.17e} {z:.17e}\n")
        for face, parent, bc in faces:
            out.write(f"{face[0]} {face[1]} {face[2]} {parent} {bc}\n")


def write_bco(path: Path) -> None:
    mappings = (506, 511, 520, 530)
    with path.open("w", encoding="utf-8") as out:
        out.write(f"{len(mappings)} 0\n")
        for bc in mappings:
            out.write(f"{bc} {bc}\n")


def write_material(path: Path, tets: Sequence[Tet]) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write(f"{len(tets)}\n")
        for ie, tet in enumerate(tets, start=1):
            out.write(f"{ie} {tet[0]} {tet[1]} {tet[2]} {tet[3]} 0\n")


def write_matprop(path: Path, viscosity: float) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("1\n")
        out.write(f"0 fluid fluid 1.0 1.0 1.0 {viscosity:.17e} 0.0\n")


def write_par(path: Path, iterations: int, re_ref: float, span: float) -> None:
    blocks = [
        ("solver_opt", "1"),
        ("restart_opt", "0"),
        ("temp_calc", "0"),
        ("freestream U V W P T", "1.0 0.0 0.0 0.0 0.0"),
        ("ntime transient_on dtfixed dtfix iwrite", f"{iterations} 0 0 1.0e-4 100"),
        ("csafm theta", "0.5 1.0"),
        (
            "CBS3D timestep controls:",
            "1 0 1 0 0 0 0 0.5 1.0e-12 1.0",
        ),
        ("re pr ra ri", f"{re_ref:.17e} 0.71 0.0 0.0"),
        ("convection_type", "0"),
        ("pnode", "1"),
        ("write_output write_time_output time_output_interval end_rtime", "0 0 1.0 1.0"),
        ("relToler absToler", "1.0e-8 1.0e-12"),
        ("velocity/temperature steady controls", "1 1.0e-7 0 1.0e-7 1 1.0e-7"),
        ("paraview tecplot nusselt_calc nusselt_flag", "1 0 0 530"),
        ("runtime_mod", "1"),
        ("alpha/source/flux controls:", "1.0 1.0 0.0 0.0"),
        ("dimensional/material/mass-flow controls:", "1 1 0"),
        ("reference/inlet/outlet controls:", f"0.0 {span:.17e} 0.0 1.0 0.0 0.0"),
        ("prescribed inlet velocity:", "1.0 0.0 0.0"),
        ("Nusselt reference controls:", "0.0 0.0 1.0"),
        ("pressure CG controls:", "1 1 5000 1.0e-10"),
        (
            "CBS3D output/monitor controls:",
            "1 10 100 0 1 500 0.0 1 20",
        ),
        ("artificial diffusion control:", "0"),
        ("Spalart-Allmaras controls:", "1 0 0"),
        ("Spalart-Allmaras inlet/thermal controls:", "3.0 0.9"),
        ("Spalart-Allmaras numerical floors:", "1.0e-14 1.0e-14 0.0"),
        ("Spalart-Allmaras switches:", "1 1 1.0e6"),
    ]

    with path.open("w", encoding="utf-8") as out:
        for label, data in blocks:
            out.write(f"{label}\n")
            out.write(f"{data}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-prefix", default="flat_plate_yplus1")
    parser.add_argument("--target-yplus", type=float, default=1.0)
    parser.add_argument("--re-x1", type=float, default=5.0e6)
    parser.add_argument("--x-min", type=float, default=-1.0 / 3.0)
    parser.add_argument("--x-max", type=float, default=2.0)
    parser.add_argument("--height", type=float, default=1.0)
    parser.add_argument("--span", type=float, default=2.0e-3)
    parser.add_argument("--nx-upstream", type=int, default=12)
    parser.add_argument("--nx-plate", type=int, default=80)
    parser.add_argument("--ny", type=int, default=48)
    parser.add_argument("--nz", type=int, default=1)
    parser.add_argument("--iterations", type=int, default=5000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.target_yplus <= 0.0 or args.re_x1 <= 0.0:
        raise SystemExit("target y+ and Reynolds number must be positive")

    nu = 1.0 / args.re_x1
    cf_ref = 0.0592 * args.re_x1 ** (-1.0 / 5.0)
    u_tau_ref = math.sqrt(0.5 * cf_ref)
    first_height = args.target_yplus * nu / u_tau_ref

    xs = streamwise_coordinates(
        args.x_min,
        args.x_max,
        args.nx_upstream,
        args.nx_plate,
    )
    ys, ratio = wall_normal_coordinates(first_height, args.height, args.ny)
    zs = [args.span * k / args.nz for k in range(args.nz + 1)]

    points, tets = make_mesh(xs, ys, zs)

    min_det = min(signed_det(points, tet) for tet in tets)
    if not (min_det > 0.0):
        raise RuntimeError(f"non-positive TET4 Jacobian detected: {min_det}")

    faces = boundary_faces(points, tets)
    counts: Dict[int, int] = {}
    for _, _, bc in faces:
        counts[bc] = counts.get(bc, 0) + 1

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    write_plt(prefix.with_suffix(".plt"), points, tets, faces)
    write_bco(prefix.with_suffix(".bco"))
    write_par(prefix.with_suffix(".par"), args.iterations, args.re_x1, args.span)
    write_material(prefix.with_suffix(".material"), tets)
    write_matprop(prefix.with_suffix(".matprop"), nu)

    metadata = {
        "reference": "NASA TMR-style zero-pressure-gradient flat plate",
        "re_x_at_x1": args.re_x1,
        "u_inf": 1.0,
        "nu_inf": nu,
        "sa_nu_tilde_inf_over_nu": 3.0,
        "target_y_plus_at_x1": args.target_yplus,
        "cf_estimate_at_x1": cf_ref,
        "u_tau_estimate_at_x1": u_tau_ref,
        "first_wall_normal_height": first_height,
        "wall_normal_growth_ratio": ratio,
        "x_min": args.x_min,
        "x_max": args.x_max,
        "height": args.height,
        "span": args.span,
        "nx_upstream": args.nx_upstream,
        "nx_plate": args.nx_plate,
        "ny": args.ny,
        "nz": args.nz,
        "nodes": len(points),
        "tetrahedra": len(tets),
        "boundary_faces": len(faces),
        "boundary_face_counts": counts,
        "minimum_signed_detJ": min_det,
    }
    prefix.with_suffix(".json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(metadata, indent=2, sort_keys=True))
    print("Run with:")
    print(f"  CBS3D_SA_WALL_TREATMENT=1 ./cbs3dpp_si {prefix}")


if __name__ == "__main__":
    main()
