#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Distributed CBS3D SA flat-plate post-processing for the TMR paper plots.

Python 3.6 compatible.  Reads a 40-rank immutable CBS partition together with
one complete distributed ASCII-VTU snapshot and writes the CSV contract used by
``tools/plot_flatplate_paper_tmr_py36.py``:

  flatplate_cf_binned.csv
  profile_x*.csv
  summary.txt

The wall-resolved skin friction is reconstructed from the P1 TET4 velocity
gradient on physical BC 530 using molecular viscosity at the wall.  Owned
points are used for profiles so MPI ghost copies are not counted twice.
"""
from __future__ import print_function

import argparse
import csv
import math
import os
import re
import xml.etree.ElementTree as ET


_NUMBER_RE = re.compile(r"[-+]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+))(?:[EeDd][-+]?\d+)?")
_STEP_RE = re.compile(r"flatplate_step_(\d{8})_rank_(\d{4})\.vtu$")


def numbers_in_line(line):
    return [m.group(0).replace("D", "E").replace("d", "e") for m in _NUMBER_RE.finditer(line)]


def read_bco(path):
    rows = []
    with open(path, "r") as handle:
        for raw in handle:
            s = raw.strip()
            if not s or s.startswith("#") or s.startswith("!"):
                continue
            rows.append(numbers_in_line(raw))
    if not rows:
        raise RuntimeError("empty bco: {}".format(path))
    nflag = int(rows[0][0])
    if len(rows) < nflag + 1:
        raise RuntimeError("truncated bco: {}".format(path))
    out = {}
    for row in rows[1:1+nflag]:
        if len(row) < 2:
            raise RuntimeError("malformed bco row in {}".format(path))
        out[int(row[0])] = int(row[1])
    return out


def read_plt(path):
    with open(path, "r") as handle:
        header = numbers_in_line(handle.readline())
        if len(header) < 3:
            raise RuntimeError("malformed plt header: {}".format(path))
        nelem, npoin, nboun = map(int, header[:3])

        tets = [None] * nelem
        for _ in range(nelem):
            row = numbers_in_line(handle.readline())
            if len(row) < 5:
                raise RuntimeError("malformed TET4 row in {}".format(path))
            ie = int(row[0])
            tets[ie-1] = tuple(map(int, row[1:5]))

        points = [None] * npoin
        for _ in range(npoin):
            row = numbers_in_line(handle.readline())
            if len(row) < 4:
                raise RuntimeError("malformed coordinate row in {}: {}".format(path, row))
            ip = int(row[0])
            points[ip-1] = (float(row[1]), float(row[2]), float(row[3]))

        faces = []
        for _ in range(nboun):
            row = numbers_in_line(handle.readline())
            if len(row) < 5:
                raise RuntimeError("malformed boundary row in {}".format(path))
            faces.append((tuple(map(int, row[:3])), int(row[3]), int(row[4])))

    if any(t is None for t in tets) or any(p is None for p in points):
        raise RuntimeError("incomplete plt indexing in {}".format(path))
    return tets, points, faces


def first_file(directory, suffix):
    names = sorted(n for n in os.listdir(directory) if n.endswith(suffix))
    if len(names) != 1:
        raise RuntimeError("expected exactly one {} in {}; found {}".format(suffix, directory, len(names)))
    return os.path.join(directory, names[0])


def read_partition(rank_dir):
    tets, points, faces = read_plt(first_file(rank_dir, ".plt"))
    mapping = read_bco(first_file(rank_dir, ".bco"))
    return tets, points, faces, mapping


def parse_data_arrays(section):
    out = {}
    if section is None:
        return out
    for array in section.findall("DataArray"):
        name = array.get("Name")
        if name and array.text:
            out[name] = [float(v) for v in array.text.split()]
    return out


def read_piece(path):
    root = ET.parse(path).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise RuntimeError("VTU Piece missing: {}".format(path))
    pnode = piece.find("./Points/DataArray")
    if pnode is None or not pnode.text:
        raise RuntimeError("VTU points missing: {}".format(path))
    raw = [float(v) for v in pnode.text.split()]
    if len(raw) % 3:
        raise RuntimeError("VTU point array malformed: {}".format(path))
    points = [(raw[i], raw[i+1], raw[i+2]) for i in range(0, len(raw), 3)]
    pd = parse_data_arrays(piece.find("PointData"))
    cd = parse_data_arrays(piece.find("CellData"))
    vel = pd.get("velocity")
    if vel is None or len(vel) != 3 * len(points):
        raise RuntimeError("velocity missing/malformed: {}".format(path))
    velocity = [(vel[i], vel[i+1], vel[i+2]) for i in range(0, len(vel), 3)]
    return points, velocity, pd, cd


def sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])


def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def tet_gradients(points, tet):
    p0, p1, p2, p3 = (points[n-1] for n in tet)
    a, b, c = sub(p1, p0), sub(p2, p0), sub(p3, p0)
    det = dot(a, cross(b, c))
    if abs(det) <= 1.0e-300:
        raise RuntimeError("singular TET4")
    g1 = tuple(x/det for x in cross(b, c))
    g2 = tuple(x/det for x in cross(c, a))
    g3 = tuple(x/det for x in cross(a, b))
    g0 = (-(g1[0]+g2[0]+g3[0]), -(g1[1]+g2[1]+g3[1]), -(g1[2]+g2[2]+g3[2]))
    return g0, g1, g2, g3


def velocity_gradient(points, tet, velocity):
    grads = tet_gradients(points, tet)
    out = []
    for comp in range(3):
        row = [0.0, 0.0, 0.0]
        for local, node in enumerate(tet):
            q = velocity[node-1][comp]
            for j in range(3):
                row[j] += q * grads[local][j]
        out.append(tuple(row))
    return tuple(out)


def complete_steps(output_dir, mpi_size):
    by_step = {}
    for name in os.listdir(output_dir):
        match = _STEP_RE.match(name)
        if not match:
            continue
        step = int(match.group(1))
        rank = int(match.group(2))
        by_step.setdefault(step, {})[rank] = os.path.join(output_dir, name)
    return dict((s, p) for s, p in by_step.items() if len(p) == mpi_size)


def write_csv(path, header, rows):
    with open(path, "w") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def interp(xs, ys, xq):
    pairs = sorted((x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y))
    if not pairs:
        return float("nan")
    if xq <= pairs[0][0]:
        return pairs[0][1]
    if xq >= pairs[-1][0]:
        return pairs[-1][1]
    for i in range(1, len(pairs)):
        x0, y0 = pairs[i-1]
        x1, y1 = pairs[i]
        if xq <= x1:
            if x1 == x0:
                return y1
            a = (xq-x0)/(x1-x0)
            return y0 + a*(y1-y0)
    return pairs[-1][1]


def encoded_x(x):
    text = ("{:.12g}".format(x)).replace("-", "m").replace(".", "p")
    return text


def main():
    parser = argparse.ArgumentParser(description="Distributed P40 flat-plate SA postprocessor")
    parser.add_argument("--partition-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--post-dir", default=None)
    parser.add_argument("--mpi-size", type=int, default=40)
    parser.add_argument("--step", type=int, default=None)
    parser.add_argument("--rho", type=float, default=1.0)
    parser.add_argument("--mu", type=float, default=2.0e-7)
    parser.add_argument("--u-inf", type=float, default=1.0)
    parser.add_argument("--wall-bc", type=int, default=530)
    parser.add_argument("--cf-bins", type=int, default=120)
    parser.add_argument("--stations", type=float, nargs="*", default=[0.97008, 1.90334])
    args = parser.parse_args()

    partition_root = os.path.abspath(args.partition_root)
    output_dir = os.path.abspath(args.output_dir)
    post_dir = os.path.abspath(args.post_dir or os.path.join(output_dir, "post_flatplate_sa_global"))
    if not os.path.isdir(partition_root):
        raise RuntimeError("partition root not found: {}".format(partition_root))
    if not os.path.isdir(output_dir):
        raise RuntimeError("output directory not found: {}".format(output_dir))
    if not os.path.isdir(post_dir):
        os.makedirs(post_dir)

    ranks = [os.path.join(partition_root, "rank_{:04d}".format(i)) for i in range(args.mpi_size)]
    missing = [p for p in ranks if not os.path.isdir(p)]
    if missing:
        raise RuntimeError("missing rank directories; first: {}".format(missing[0]))

    steps = complete_steps(output_dir, args.mpi_size)
    if not steps:
        raise RuntimeError("no complete {}-piece VTU snapshot under {}".format(args.mpi_size, output_dir))
    step = args.step if args.step is not None else max(steps)
    if step not in steps:
        raise RuntimeError("requested step {} is not a complete distributed snapshot".format(step))

    wall_faces = []
    seen_faces = set()
    owned_nodes = {}
    max_abs_w = 0.0
    max_mut = 0.0

    for irank, rank_dir in enumerate(ranks):
        tets, mesh_points, boundary, mapping = read_partition(rank_dir)
        vtu_points, velocity, pd, cd = read_piece(steps[step][irank])
        if len(mesh_points) != len(vtu_points):
            raise RuntimeError("mesh/VTU point count mismatch rank {}".format(irank))

        gids = [int(round(v)) for v in pd.get("global_node_id", [])]
        owned = [int(round(v)) for v in pd.get("is_owned", [])]
        nu_tilde = pd.get("nu_tilde")
        mut = pd.get("mu_t_over_mu")
        wall_distance = pd.get("wall_distance")
        if not (len(gids) == len(vtu_points) == len(owned)):
            raise RuntimeError("ownership/global-node fields missing on rank {}".format(irank))
        if nu_tilde is None or mut is None or wall_distance is None:
            raise RuntimeError("SA point fields missing on rank {}".format(irank))

        for i, gid in enumerate(gids):
            if owned[i] == 1:
                owned_nodes[gid] = (vtu_points[i], velocity[i], nu_tilde[i], mut[i], wall_distance[i])
                max_abs_w = max(max_abs_w, abs(velocity[i][2]))
                max_mut = max(max_mut, mut[i])

        for face, parent, raw_bc in boundary:
            if mapping.get(raw_bc) != args.wall_bc:
                continue
            pts = [mesh_points[n-1] for n in face]
            xc = sum(p[0] for p in pts) / 3.0
            if xc <= 0.0:
                continue
            gface = tuple(sorted(gids[n-1] for n in face))
            if gface in seen_faces:
                continue
            seen_faces.add(gface)

            tet = tets[parent-1]
            opp = [n for n in tet if n not in set(face)]
            if len(opp) != 1:
                raise RuntimeError("wall face does not match parent TET4")
            opposite = mesh_points[opp[0]-1]
            raw_n = cross(sub(pts[1], pts[0]), sub(pts[2], pts[0]))
            nmag = norm(raw_n)
            if not nmag > 0.0:
                raise RuntimeError("degenerate wall face")
            area = 0.5 * nmag
            n = tuple(v/nmag for v in raw_n)
            centroid = tuple(sum(p[j] for p in pts)/3.0 for j in range(3))
            if dot(n, sub(opposite, centroid)) > 0.0:
                n = (-n[0], -n[1], -n[2])

            grad = velocity_gradient(mesh_points, tet, velocity)
            traction = tuple((args.mu/args.rho) * dot(grad[i], n) for i in range(3))
            normal_component = dot(traction, n)
            tangential = tuple(traction[i] - normal_component*n[i] for i in range(3))
            shear_over_rho = norm(tangential)
            cf = 2.0 * shear_over_rho / (args.u_inf * args.u_inf)
            tau_w = args.rho * shear_over_rho
            wall_faces.append((xc, area, cf, tau_w))

    if not wall_faces:
        raise RuntimeError("no physical BC {} wall faces found".format(args.wall_bc))
    if not owned_nodes:
        raise RuntimeError("no owned distributed nodes found")

    xmin = min(r[0] for r in wall_faces)
    xmax = max(r[0] for r in wall_faces)
    edges = [xmin + (xmax-xmin)*i/float(args.cf_bins) for i in range(args.cf_bins+1)]
    cf_rows = []
    for i in range(args.cf_bins):
        lo, hi = edges[i], edges[i+1]
        subset = [r for r in wall_faces if (r[0] >= lo and (r[0] < hi or (i == args.cf_bins-1 and r[0] <= hi)))]
        if not subset:
            continue
        area = sum(r[1] for r in subset)
        xbar = sum(r[0]*r[1] for r in subset) / area
        cf = sum(r[2]*r[1] for r in subset) / area
        tau = sum(r[3]*r[1] for r in subset) / area
        cf_rows.append((xbar, cf, tau, len(subset), area))

    cf_path = os.path.join(post_dir, "flatplate_cf_binned.csv")
    write_csv(cf_path, ["x", "Cf", "tau_w", "face_count", "area"], cf_rows)

    cf_x = [r[0] for r in cf_rows]
    cf_y = [r[1] for r in cf_rows]
    nu = args.mu / args.rho
    profile_paths = []

    unique_x = sorted(set(round(v[0][0], 12) for v in owned_nodes.values() if v[0][0] >= 0.0))
    for target in args.stations:
        actual_x = min(unique_x, key=lambda x: abs(x-target))
        selected = [v for v in owned_nodes.values() if abs(v[0][0]-actual_x) <= 2.0e-11]
        by_y = {}
        for coord, vel, q, ratio, wd in selected:
            by_y.setdefault(round(coord[1], 14), []).append((coord, vel, q, ratio, wd))

        local_cf = interp(cf_x, cf_y, actual_x)
        if not local_cf > 0.0:
            raise RuntimeError("non-positive interpolated Cf at x={}".format(actual_x))
        u_tau = args.u_inf * math.sqrt(0.5 * local_cf)

        rows = []
        for ykey in sorted(by_y):
            vals = by_y[ykey]
            n = float(len(vals))
            x = sum(v[0][0] for v in vals)/n
            y = sum(v[0][1] for v in vals)/n
            z = sum(v[0][2] for v in vals)/n
            u = sum(v[1][0] for v in vals)/n
            vv = sum(v[1][1] for v in vals)/n
            w = sum(v[1][2] for v in vals)/n
            vm = sum(norm(v[1]) for v in vals)/n
            q = sum(v[2] for v in vals)/n
            ratio = sum(v[3] for v in vals)/n
            wd = sum(v[4] for v in vals)/n
            yp = wd * u_tau / nu if wd >= 0.0 else float("nan")
            up = u / u_tau
            rows.append((x, y, z, u, vv, w, vm, q, ratio, wd, yp, up))

        name = "profile_x{}.csv".format(encoded_x(actual_x))
        path = os.path.join(post_dir, name)
        write_csv(path,
                  ["x", "y", "z", "u", "v", "w", "velocity_magnitude", "nu_tilde", "mu_t_over_mu", "wall_distance", "y_plus", "u_plus"],
                  rows)
        profile_paths.append(path)

    with open(os.path.join(post_dir, "summary.txt"), "w") as handle:
        handle.write("iteration = {}\n".format(step))
        handle.write("mpi_size = {}\n".format(args.mpi_size))
        handle.write("wall_faces_unique = {}\n".format(len(wall_faces)))
        handle.write("owned_global_nodes = {}\n".format(len(owned_nodes)))
        handle.write("max_abs_w_over_all_piece_nodes = {:.17g}\n".format(max_abs_w))
        handle.write("max_mu_t_over_mu_owned_nodes = {:.17g}\n".format(max_mut))
        handle.write("rho = {:.17g}\n".format(args.rho))
        handle.write("mu = {:.17g}\n".format(args.mu))
        handle.write("u_inf = {:.17g}\n".format(args.u_inf))

    print("="*72)
    print("DISTRIBUTED FLAT-PLATE CSV POSTPROCESS: PASS")
    print("="*72)
    print("step       : {}".format(step))
    print("wall faces : {}".format(len(wall_faces)))
    print("owned nodes: {}".format(len(owned_nodes)))
    print("Cf CSV     : {}".format(cf_path))
    for path in profile_paths:
        print("profile    : {}".format(path))
    print("summary    : {}".format(os.path.join(post_dir, "summary.txt")))
    print("post-dir   : {}".format(post_dir))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("ERROR: {}".format(exc))
        raise
