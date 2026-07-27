#!/usr/bin/env python3
"""Audit explicit CBS momentum/energy timestep limits from rank-local VTU output.

The script is read-only. It reproduces the geometric timestep estimate currently
used by TimeStep.cpp and also assembles conservative lumped-capacity thermal
forward-Euler bounds from the actual P1 tetrahedral mesh.
"""

import argparse
import math
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


def fail(message):
    raise RuntimeError(message)


def floats(node, context):
    try:
        return [float(v) for v in (node.text or "").split()]
    except ValueError as exc:
        raise RuntimeError("{}: invalid floating-point DataArray".format(context)) from exc


def ints(node, context):
    try:
        return [int(v) for v in (node.text or "").split()]
    except ValueError as exc:
        raise RuntimeError("{}: invalid integer DataArray".format(context)) from exc


def named(parent, candidates, context):
    wanted = set(candidates)
    for node in parent.findall("DataArray"):
        if node.get("Name") in wanted:
            return node
    fail("{}: missing DataArray; expected one of {}".format(context, sorted(wanted)))


def vec_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a):
    return math.sqrt(dot(a, a))


def triangle_area(a, b, c):
    return 0.5 * norm(cross(vec_sub(b, a), vec_sub(c, a)))


def gradients_and_volume(p0, p1, p2, p3):
    a = vec_sub(p1, p0)
    b = vec_sub(p2, p0)
    c = vec_sub(p3, p0)
    det = dot(a, cross(b, c))
    if det == 0.0 or not math.isfinite(det):
        fail("degenerate tetrahedron encountered")

    # J has columns a,b,c. Rows of J^{-1} are grad(N1), grad(N2), grad(N3).
    g1 = tuple(v / det for v in cross(b, c))
    g2 = tuple(v / det for v in cross(c, a))
    g3 = tuple(v / det for v in cross(a, b))
    g0 = (-(g1[0] + g2[0] + g3[0]),
          -(g1[1] + g2[1] + g3[1]),
          -(g1[2] + g2[2] + g3[2]))
    return (g0, g1, g2, g3), abs(det) / 6.0


def read_matprop(path):
    tokens = path.read_text(encoding="utf-8").split()
    if not tokens:
        fail("empty material-property file: {}".format(path))
    nmat = int(tokens[0])
    expected = 1 + 8 * nmat
    if len(tokens) < expected:
        fail("{}: incomplete material-property records".format(path))

    materials = {}
    pos = 1
    for _ in range(nmat):
        mid = int(tokens[pos])
        rho = float(tokens[pos + 3])
        cp = float(tokens[pos + 4])
        k = float(tokens[pos + 5])
        mu = float(tokens[pos + 6])
        pos += 8
        if rho <= 0.0 or cp <= 0.0 or k <= 0.0 or mu < 0.0:
            fail("non-physical material properties for ID {}".format(mid))
        materials[mid] = {
            "rho": rho,
            "cp": cp,
            "rho_cp": rho * cp,
            "k": k,
            "mu": mu,
            "alpha": k / (rho * cp),
            "nu": mu / rho,
        }
    return materials


def read_piece(path):
    root = ET.parse(path).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        fail("{}: missing UnstructuredGrid/Piece".format(path))

    points_parent = piece.find("Points")
    cells_parent = piece.find("Cells")
    point_data = piece.find("PointData")
    cell_data = piece.find("CellData")
    if any(node is None for node in (points_parent, cells_parent, point_data, cell_data)):
        fail("{}: incomplete VTU Piece".format(path))

    flat_points = floats(points_parent.find("DataArray"), str(path))
    if len(flat_points) % 3:
        fail("{}: invalid point coordinate length".format(path))
    points = [tuple(flat_points[i:i + 3]) for i in range(0, len(flat_points), 3)]

    flat_velocity = floats(named(point_data, ("velocity",), str(path)), str(path))
    if len(flat_velocity) != 3 * len(points):
        fail("{}: velocity array length mismatch".format(path))
    velocity = [tuple(flat_velocity[i:i + 3]) for i in range(0, len(flat_velocity), 3)]

    global_node_id = ints(named(point_data, ("global_node_id",), str(path)), str(path))
    connectivity = ints(named(cells_parent, ("connectivity",), str(path)), str(path))
    offsets = ints(named(cells_parent, ("offsets",), str(path)), str(path))
    types = ints(named(cells_parent, ("types",), str(path)), str(path))
    material_id = ints(named(cell_data, ("material_id", "mat_id"), str(path)), str(path))
    global_element_id = ints(
        named(cell_data, ("global_element_id", "parent_element"), str(path)), str(path)
    )

    ncells = len(offsets)
    if not (len(types) == len(material_id) == len(global_element_id) == ncells):
        fail("{}: cell array length mismatch".format(path))
    if len(global_node_id) != len(points):
        fail("{}: global_node_id length mismatch".format(path))

    return points, velocity, global_node_id, connectivity, offsets, types, material_id, global_element_id


def update_min(record, value, details):
    if value < record[0]:
        record[0] = value
        record[1] = details


def fmt_control(details):
    if details is None:
        return "n/a"
    return (
        "global_element={gid} material={mid} h={h:.12e} m "
        "U={speed:.12e} m/s D={diff:.12e} m^2/s file={file}"
    ).format(**details)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtu-dir", type=Path, required=True)
    parser.add_argument("--step", type=int, required=True)
    parser.add_argument("--matprop", type=Path, required=True)
    parser.add_argument("--csafm", type=float, default=0.70)
    parser.add_argument("--current-dt", type=float, default=1.0e-5)
    args = parser.parse_args()

    if args.step < 0 or args.csafm <= 0.0 or args.current_dt <= 0.0:
        fail("step, csafm and current-dt must be positive")

    materials = read_matprop(args.matprop)
    tag = "{:08d}".format(args.step)
    files = sorted(p for p in args.vtu_dir.glob("*.vtu") if tag in p.name)
    if not files:
        fail("no VTU files found for step {}".format(args.step))

    capacity = defaultdict(float)
    kdiag = defaultdict(float)
    krow_abs_element_bound = defaultdict(float)

    min_h = [math.inf, None]
    min_h_fluid = [math.inf, None]
    min_h_solid = [math.inf, None]
    min_adv = [math.inf, None]
    min_diff = [math.inf, None]
    min_diff_fluid = [math.inf, None]
    min_diff_solid = [math.inf, None]
    min_code = [math.inf, None]

    tetra_count = 0
    volume_total = 0.0

    for piece_index, path in enumerate(files, 1):
        print("Reading {:2d}/{:2d}: {}".format(piece_index, len(files), path.name), flush=True)
        (points, velocity, gids, conn, offsets, types, mids, geids) = read_piece(path)
        start = 0
        for cell_index, end in enumerate(offsets):
            nodes = conn[start:end]
            start = end
            if types[cell_index] != 10:
                continue
            if len(nodes) != 4:
                fail("{}: tetrahedron with {} nodes".format(path, len(nodes)))

            mid = mids[cell_index]
            if mid not in materials:
                fail("{}: material ID {} missing from matprop".format(path, mid))
            mat = materials[mid]

            local_points = [points[i] for i in nodes]
            grads, volume = gradients_and_volume(*local_points)
            if volume <= 0.0 or not math.isfinite(volume):
                fail("{}: invalid tetrahedral volume".format(path))

            face_areas = (
                triangle_area(local_points[1], local_points[2], local_points[3]),
                triangle_area(local_points[0], local_points[3], local_points[2]),
                triangle_area(local_points[0], local_points[1], local_points[3]),
                triangle_area(local_points[0], local_points[2], local_points[1]),
            )
            if min(face_areas) <= 0.0:
                fail("{}: invalid tetrahedral face area".format(path))
            h = min(3.0 * volume / area for area in face_areas)

            speed = max(norm(velocity[i]) for i in nodes)
            diff = mat["alpha"] if mid != 0 else max(mat["alpha"], mat["nu"])
            dt_adv = math.inf if mid != 0 else h / max(speed, 1.0e-14)
            dt_diff = h * h / (2.0 * diff)
            dt_code = args.csafm * min(dt_adv, dt_diff)

            details = {
                "gid": geids[cell_index],
                "mid": mid,
                "h": h,
                "speed": speed,
                "diff": diff,
                "file": path.name,
            }
            update_min(min_h, h, details)
            update_min(min_code, dt_code, details)
            update_min(min_diff, dt_diff, details)
            if mid == 0:
                update_min(min_h_fluid, h, details)
                update_min(min_adv, dt_adv, details)
                update_min(min_diff_fluid, dt_diff, details)
            else:
                update_min(min_h_solid, h, details)
                update_min(min_diff_solid, dt_diff, details)

            cap_local = mat["rho_cp"] * volume / 4.0
            kval = mat["k"] * volume
            for a in range(4):
                gid = gids[nodes[a]]
                capacity[gid] += cap_local
                row_abs = 0.0
                for b in range(4):
                    kab = kval * dot(grads[a], grads[b])
                    row_abs += abs(kab)
                    if a == b:
                        kdiag[gid] += kab
                # Sum absolute element rows. This is conservative because it
                # intentionally ignores cancellation during global assembly.
                krow_abs_element_bound[gid] += row_abs

            tetra_count += 1
            volume_total += volume

    lambda_upper = 0.0
    dt_self = math.inf
    controlling_lambda_gid = None
    controlling_self_gid = None
    for gid, cap in capacity.items():
        if cap <= 0.0:
            fail("non-positive assembled thermal capacity at global node {}".format(gid))
        row = krow_abs_element_bound[gid]
        diag = kdiag[gid]
        lam = row / cap
        if lam > lambda_upper:
            lambda_upper = lam
            controlling_lambda_gid = gid
        if diag > 0.0:
            candidate = cap / diag
            if candidate < dt_self:
                dt_self = candidate
                controlling_self_gid = gid

    dt_fe_conservative = 2.0 / lambda_upper if lambda_upper > 0.0 else math.inf

    print("\n=== Actual-mesh CBS timestep audit ===")
    print("VTU pieces                         : {}".format(len(files)))
    print("Owned tetrahedra processed         : {}".format(tetra_count))
    print("Assembled unique thermal nodes     : {}".format(len(capacity)))
    print("Total tetrahedral volume           : {:.12e} m^3".format(volume_total))
    print("CSAFM                               : {:.6f}".format(args.csafm))
    print("Current fixed dt                    : {:.12e} s".format(args.current_dt))

    print("\nMinimum tetrahedral altitude        : {:.12e} m".format(min_h[0]))
    print("  controller: {}".format(fmt_control(min_h[1])))
    print("Minimum fluid altitude              : {:.12e} m".format(min_h_fluid[0]))
    print("Minimum solid altitude              : {:.12e} m".format(min_h_solid[0]))

    print("\nMinimum raw fluid advective h/U     : {:.12e} s".format(min_adv[0]))
    print("  controller: {}".format(fmt_control(min_adv[1])))
    print("Minimum raw diffusion h^2/(2D)      : {:.12e} s".format(min_diff[0]))
    print("  controller: {}".format(fmt_control(min_diff[1])))
    print("Minimum raw fluid diffusion limit   : {:.12e} s".format(min_diff_fluid[0]))
    print("Minimum raw solid diffusion limit   : {:.12e} s".format(min_diff_solid[0]))

    print("\nCode-equivalent global dt estimate  : {:.12e} s".format(min_code[0]))
    print("  controller: {}".format(fmt_control(min_code[1])))
    print("  margin over current dt            : {:.6f} x".format(min_code[0] / args.current_dt))

    print("\nThermal FE conservative Euler bound : {:.12e} s".format(dt_fe_conservative))
    print("  derived from 2/max row-bound      : global_node={}".format(controlling_lambda_gid))
    print("  margin over current dt            : {:.6f} x".format(dt_fe_conservative / args.current_dt))
    print("Thermal nonnegative self-weight dt  : {:.12e} s".format(dt_self))
    print("  controller                        : global_node={}".format(controlling_self_gid))
    print("  margin over current dt            : {:.6f} x".format(dt_self / args.current_dt))

    print("\nInterpretation:")
    print("  - The code-equivalent value reproduces the current h/U and h^2/(2D) heuristic.")
    print("  - The FE bound uses the assembled lumped capacity and conductivity on the actual mesh.")
    print("  - Neither bound includes nonlinear transient accuracy requirements; this is a steady pseudo-time audit.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print("ERROR: {}".format(exc), file=__import__("sys").stderr)
        raise SystemExit(1)
