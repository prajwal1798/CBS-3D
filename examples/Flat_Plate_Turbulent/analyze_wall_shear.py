#!/usr/bin/env python3
"""Extract flat-plate wall shear, Cf, u_tau and y+ from CBS3D output.

The script uses only the Python standard library.  It reads the solver-native
.plt/.bco files and one ASCII VTU solution.  Two shear definitions are supported:

  wall-model : reproduce the production Spalding traction from the parent TET4
               opposite-node velocity and tetrahedron altitude.
  resolved   : evaluate the molecular viscous natural traction nu grad(u).n from
               the constant P1 TET4 velocity gradient.

The resulting face CSV is suitable for comparison with NASA TMR flat-plate SA
surface skin-friction data.  Correlation columns are included only as secondary
engineering checks, not as substitutes for TMR verification data.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

Vec = Tuple[float, float, float]
Tet = Tuple[int, int, int, int]
Face = Tuple[int, int, int]


def add(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def mul(a: Vec, scalar: float) -> Vec:
    return (a[0] * scalar, a[1] * scalar, a[2] * scalar)


def dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Vec, b: Vec) -> Vec:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a: Vec) -> float:
    return math.sqrt(dot(a, a))


def unit(a: Vec) -> Vec:
    magnitude = norm(a)
    if not magnitude > 0.0:
        raise ValueError("zero vector cannot be normalized")
    return mul(a, 1.0 / magnitude)


def tangent(v: Vec, normal: Vec) -> Vec:
    n = unit(normal)
    return sub(v, mul(n, dot(v, n)))


def triangle_geometry(points: Sequence[Vec], face: Face, opposite: Vec) -> Tuple[Vec, float, Vec]:
    p0, p1, p2 = (points[node - 1] for node in face)
    raw = cross(sub(p1, p0), sub(p2, p0))
    raw_norm = norm(raw)
    if not raw_norm > 0.0:
        raise ValueError("degenerate wall triangle")
    area = 0.5 * raw_norm
    normal = mul(raw, 1.0 / raw_norm)
    centroid = mul(add(add(p0, p1), p2), 1.0 / 3.0)

    # Outward normal points away from the parent tetrahedron's opposite vertex.
    if dot(normal, sub(opposite, centroid)) > 0.0:
        normal = mul(normal, -1.0)
    return centroid, area, normal


def det3(a: Vec, b: Vec, c: Vec) -> float:
    return dot(a, cross(b, c))


def inverse3(columns: Tuple[Vec, Vec, Vec]) -> Tuple[Vec, Vec, Vec]:
    a, b, c = columns
    determinant = det3(a, b, c)
    if abs(determinant) <= 1.0e-300:
        raise ValueError("singular tetrahedron Jacobian")

    # Rows of J^-1 for J=[a b c].
    r0 = mul(cross(b, c), 1.0 / determinant)
    r1 = mul(cross(c, a), 1.0 / determinant)
    r2 = mul(cross(a, b), 1.0 / determinant)
    return r0, r1, r2


def tet_gradients(points: Sequence[Vec], tet: Tet) -> Tuple[Vec, Vec, Vec, Vec]:
    p0, p1, p2, p3 = (points[node - 1] for node in tet)
    inv_rows = inverse3((sub(p1, p0), sub(p2, p0), sub(p3, p0)))
    g1, g2, g3 = inv_rows
    g0 = mul(add(add(g1, g2), g3), -1.0)
    return g0, g1, g2, g3


def velocity_gradient(points: Sequence[Vec], tet: Tet, velocity: Sequence[Vec]) -> Tuple[Vec, Vec, Vec]:
    gradients = tet_gradients(points, tet)
    rows: List[Vec] = []
    for component in range(3):
        row = [0.0, 0.0, 0.0]
        for local, node in enumerate(tet):
            value = velocity[node - 1][component]
            for dim in range(3):
                row[dim] += value * gradients[local][dim]
        rows.append((row[0], row[1], row[2]))
    return rows[0], rows[1], rows[2]


def read_bco(path: Path) -> Dict[int, int]:
    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith(("#", "!"))
    ]
    if not lines:
        raise ValueError("empty .bco file")
    nflag = int(lines[0].split()[0])
    if len(lines) < nflag + 1:
        raise ValueError("truncated .bco file")
    mapping: Dict[int, int] = {}
    for line in lines[1 : nflag + 1]:
        fields = line.split()
        mapping[int(fields[0])] = int(fields[1])
    return mapping


def read_plt(path: Path) -> Tuple[List[Tet], List[Vec], List[Tuple[Face, int, int]]]:
    tokens = path.read_text(encoding="utf-8").split()
    cursor = 0

    def take_int() -> int:
        nonlocal cursor
        value = int(tokens[cursor])
        cursor += 1
        return value

    def take_float() -> float:
        nonlocal cursor
        value = float(tokens[cursor])
        cursor += 1
        return value

    nelem, npoin, nboun = take_int(), take_int(), take_int()
    tets: List[Tet] = [(-1, -1, -1, -1)] * nelem
    for _ in range(nelem):
        ie = take_int()
        tet = (take_int(), take_int(), take_int(), take_int())
        tets[ie - 1] = tet

    points: List[Vec] = [(0.0, 0.0, 0.0)] * npoin
    for _ in range(npoin):
        ip = take_int()
        points[ip - 1] = (take_float(), take_float(), take_float())

    faces: List[Tuple[Face, int, int]] = []
    for _ in range(nboun):
        face = (take_int(), take_int(), take_int())
        parent = take_int()
        raw_bc = take_int()
        faces.append((face, parent, raw_bc))

    if cursor != len(tokens):
        raise ValueError("unexpected trailing tokens in .plt file")
    return tets, points, faces


def read_vtu_velocity(path: Path) -> Tuple[List[Vec], List[Vec]]:
    root = ET.parse(path).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError("VTU Piece not found")

    points_array = piece.find("./Points/DataArray")
    if points_array is None or not points_array.text:
        raise ValueError("VTU point coordinates not found")
    point_values = [float(v) for v in points_array.text.split()]
    if len(point_values) % 3:
        raise ValueError("invalid VTU point coordinate count")
    points = [
        (point_values[i], point_values[i + 1], point_values[i + 2])
        for i in range(0, len(point_values), 3)
    ]

    velocity_array = None
    for array in piece.findall("./PointData/DataArray"):
        if array.get("Name") == "velocity":
            velocity_array = array
            break
    if velocity_array is None or not velocity_array.text:
        raise ValueError("VTU velocity field not found")
    values = [float(v) for v in velocity_array.text.split()]
    if len(values) != 3 * len(points):
        raise ValueError("VTU velocity size does not match point count")
    velocity = [
        (values[i], values[i + 1], values[i + 2])
        for i in range(0, len(values), 3)
    ]
    return points, velocity


def exponential_remainder(x: float, order: int) -> float:
    if x < 1.0:
        first = order + 1
        term = 1.0
        for n in range(1, first + 1):
            term *= x / n
        total = term
        for n in range(first + 1, 101):
            term *= x / n
            total += term
            if abs(term) <= 1.0e-16 * max(1.0, abs(total)):
                break
        return total
    polynomial = 1.0 + x + 0.5 * x * x + x**3 / 6.0
    if order == 4:
        polynomial += x**4 / 24.0
    return math.exp(x) - polynomial


def spalding_yplus(uplus: float, kappa: float, b: float) -> float:
    x = kappa * uplus
    return uplus + math.exp(-kappa * b) * exponential_remainder(x, 4)


def solve_spalding(speed: float, wall_height: float, nu: float, kappa: float, b: float) -> Tuple[float, float, float]:
    if speed <= 1.0e-14:
        return 0.0, 0.0, 0.0
    rey = speed * wall_height / nu
    lo, hi = 0.0, math.sqrt(rey)
    uplus = math.log(max(1.0, rey)) / kappa + b
    if not lo < uplus < hi:
        uplus = 0.5 * (lo + hi)

    coefficient = math.exp(-kappa * b)
    for _ in range(100):
        x = kappa * uplus
        yplus = uplus + coefficient * exponential_remainder(x, 4)
        derivative_y = 1.0 + coefficient * kappa * exponential_remainder(x, 3)
        residual = uplus * yplus - rey
        if abs(residual) <= 1.0e-13 * rey:
            utau = speed / uplus
            return utau, uplus, yplus
        if residual < 0.0:
            lo = uplus
        else:
            hi = uplus
        derivative = yplus + uplus * derivative_y
        candidate = uplus - residual / derivative
        if not math.isfinite(candidate) or not lo < candidate < hi:
            candidate = 0.5 * (lo + hi)
        uplus = candidate
    raise RuntimeError("Spalding root solve did not converge")


def find_opposite(tet: Tet, face: Face) -> int:
    face_set = set(face)
    candidates = [node for node in tet if node not in face_set]
    if len(candidates) != 1:
        raise ValueError("boundary face is not a face of its parent tetrahedron")
    return candidates[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_prefix", type=Path)
    parser.add_argument("vtu", type=Path)
    parser.add_argument("--mode", choices=("wall-model", "resolved"), default="wall-model")
    parser.add_argument("--wall-bc", type=int, default=530)
    parser.add_argument("--re-x1", type=float, default=None)
    parser.add_argument("--u-inf", type=float, default=1.0)
    parser.add_argument("--kappa", type=float, default=0.41)
    parser.add_argument("--B", type=float, default=5.2)
    parser.add_argument("--csv", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    prefix = args.case_prefix
    plt_path = prefix.with_suffix(".plt")
    bco_path = prefix.with_suffix(".bco")
    manifest_path = prefix.with_suffix(".json")

    re_x1 = args.re_x1
    if re_x1 is None and manifest_path.exists():
        re_x1 = float(json.loads(manifest_path.read_text(encoding="utf-8"))["re_x_at_x1"])
    if re_x1 is None or re_x1 <= 0.0:
        raise SystemExit("Re_x(x=1) must be supplied by --re-x1 or the generator JSON manifest")
    if args.u_inf <= 0.0:
        raise SystemExit("U_inf must be positive")
    nu = args.u_inf / re_x1

    tets, mesh_points, boundary = read_plt(plt_path)
    mapping = read_bco(bco_path)
    vtu_points, velocity = read_vtu_velocity(args.vtu)
    if len(vtu_points) != len(mesh_points):
        raise RuntimeError("VTU and .plt node counts differ")

    coordinate_error = max(
        norm(sub(a, b)) for a, b in zip(mesh_points, vtu_points)
    )
    if coordinate_error > 1.0e-10:
        raise RuntimeError(f"VTU/.plt coordinates disagree: max error={coordinate_error:.6e}")

    rows: List[Dict[str, float]] = []
    total_drag = 0.0
    total_area = 0.0

    for face, parent, raw_bc in boundary:
        if mapping.get(raw_bc) != args.wall_bc:
            continue
        tet = tets[parent - 1]
        opposite_node = find_opposite(tet, face)
        opposite = mesh_points[opposite_node - 1]
        centroid, area, normal = triangle_geometry(mesh_points, face, opposite)

        if args.mode == "wall-model":
            sample = velocity[opposite_node - 1]
            sample_tangent = tangent(sample, normal)
            sample_speed = norm(sample_tangent)
            # TET4 altitude h = 3V/A = |detJ|/(2A).
            p0, p1, p2, p3 = (mesh_points[node - 1] for node in tet)
            detj = abs(det3(sub(p1, p0), sub(p2, p0), sub(p3, p0)))
            sample_height = detj / (2.0 * area)
            utau, uplus, yplus = solve_spalding(
                sample_speed,
                sample_height,
                nu,
                args.kappa,
                args.B,
            )
            shear_mag = utau * utau
        else:
            grad_u = velocity_gradient(mesh_points, tet, velocity)
            traction = (
                nu * dot(grad_u[0], normal),
                nu * dot(grad_u[1], normal),
                nu * dot(grad_u[2], normal),
            )
            tangential_traction = tangent(traction, normal)
            shear_mag = norm(tangential_traction)
            utau = math.sqrt(max(0.0, shear_mag))
            sample_height = 0.0
            uplus = 0.0
            yplus = 0.0
            if utau > 0.0:
                # Report the opposite-node height in wall units even for the
                # resolved mode, for direct mesh-resolution auditing.
                p0, p1, p2, p3 = (mesh_points[node - 1] for node in tet)
                detj = abs(det3(sub(p1, p0), sub(p2, p0), sub(p3, p0)))
                sample_height = detj / (2.0 * area)
                yplus = sample_height * utau / nu
                sample_speed = norm(tangent(velocity[opposite_node - 1], normal))
                uplus = sample_speed / utau

        cf = 2.0 * shear_mag / (args.u_inf * args.u_inf)
        rex = re_x1 * max(centroid[0], 0.0)
        cf_one_fifth = 0.0592 * rex ** (-1.0 / 5.0) if rex > 0.0 else math.nan
        cf_one_seventh = 0.025 * rex ** (-1.0 / 7.0) if rex > 0.0 else math.nan

        rows.append(
            {
                "x": centroid[0],
                "y": centroid[1],
                "z": centroid[2],
                "area": area,
                "parent_element": float(parent),
                "opposite_node": float(opposite_node),
                "sample_height": sample_height,
                "u_tau": utau,
                "u_plus_sample": uplus,
                "y_plus_sample": yplus,
                "cf": cf,
                "cf_power_1_5": cf_one_fifth,
                "cf_power_1_7": cf_one_seventh,
            }
        )
        total_drag += shear_mag * area
        total_area += area

    if not rows:
        raise RuntimeError(f"no boundary faces mapped to wall BC {args.wall_bc}")

    rows.sort(key=lambda row: (row["x"], row["z"]))
    csv_path = args.csv or prefix.with_name(prefix.name + f"_{args.mode}_wall.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    positive_x = [row for row in rows if row["x"] > 0.0]
    target = min(positive_x, key=lambda row: abs(row["x"] - 0.97008))
    min_yplus = min(row["y_plus_sample"] for row in positive_x)
    max_yplus = max(row["y_plus_sample"] for row in positive_x)
    area_weighted_cf = (
        sum(row["cf"] * row["area"] for row in positive_x)
        / sum(row["area"] for row in positive_x)
    )

    summary = {
        "mode": args.mode,
        "wall_faces": len(rows),
        "wall_area": total_area,
        "kinematic_drag_integral": total_drag,
        "area_weighted_cf": area_weighted_cf,
        "min_sample_y_plus": min_yplus,
        "max_sample_y_plus": max_yplus,
        "nearest_x_to_0.97008": target["x"],
        "cf_near_x_0.97008": target["cf"],
        "y_plus_near_x_0.97008": target["y_plus_sample"],
        "u_tau_near_x_0.97008": target["u_tau"],
        "csv": str(csv_path),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
