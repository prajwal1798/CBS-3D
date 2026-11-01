#!/usr/bin/env python3
# ================================
# FILE      : solid_energy.py
# PURPOSE   : SOLID THERMAL ENERGY OF A CBS3D CHT VTU, FOR THE TRANSIENT
#             ENERGY BUDGET OF A TIME-AVERAGED RUN.
#
#             OVER AN AVERAGING WINDOW OF LENGTH dt THE GLOBAL BALANCE IS
#
#                 <Q_out> = Q_in - dE_solid/dt - dE_fluid/dt
#
#             WITH dE_fluid/dt NEGLIGIBLE. THE PRODUCT-OF-MEANS TOOL
#             energy_balance.py --mean EVALUATES <u><T>, SO ITS DEFICIT
#             CONTAINS BOTH THE SOLID STORAGE AND THE <u'T'> COVARIANCE.
#             THIS TOOL MEASURES THE STORAGE TERM FROM SAVED SNAPSHOTS.
#
#             EVERY SOLID TETRAHEDRON IS COUNTED ONCE THROUGH
#             global_element_id (THE DISTRIBUTED VTU REPEATS GHOST CELLS).
#             E_s = rho cp SUM_e V_e (T1+T2+T3+T4)/4   (EXACT FOR P1 T)
#
# USAGE     : python3 tools/solid_energy.py RESULT [RESULT ...]
#                 [--rho 7850] [--cp 500] [--solid-id 1]
#             RESULT = .pvtu / .vtu / FOLDER (LATEST FILE UNDER IT)
#
# OUTPUT    : ONE LINE PER SNAPSHOT: FILE, SOLID CELLS, VOLUME, MEAN T, E_s.
#             SUBTRACT TWO E_s VALUES AND DIVIDE BY THE TIME INTERVAL TO GET
#             THE STORAGE RATE. EXIT 0 ON SUCCESS.
#
# STATUS    : NUMPY ONLY. CBS3D ASCII VTU WITH CellData material_id AND
#             global_element_id.
# ================================

from __future__ import print_function

import argparse
import os
import re
import sys

import numpy as np

TET = 10


def _array(text, name):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                  text, re.S)
    if m is None:
        raise ValueError("NO DataArray '%s'" % name)
    return np.array(m.group(1).split(), dtype=float)


def _array_cell(text, name):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                  text, re.S)
    if m is None:
        return None
    return np.array(m.group(1).split(), dtype=np.int64)


def latest(path):
    if not os.path.isdir(path):
        return path
    found = []
    for root, _dirs, files in os.walk(path):
        for f in files:
            if f.endswith(".pvtu") or (f.endswith(".vtu") and "_rank_" not in f
                                       and "boundary" not in f):
                full = os.path.join(root, f)
                found.append((os.path.getmtime(full), full))
    if not found:
        raise ValueError("NO .pvtu OR .vtu UNDER %s" % path)
    return max(found)[1]


def read_piece(path):
    """POINTS, SOLID TETS, T AT THOSE NODES, GLOBAL CELL IDS."""
    with open(path) as handle:
        text = handle.read()
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    points = np.array(m.group(1).split(), dtype=float).reshape(-1, 3)
    conn = _array(text, "connectivity").astype(np.int64)
    offsets = _array(text, "offsets").astype(np.int64)
    types = _array(text, "types").astype(np.int64)
    temperature = _array(text, "temperature")
    material = _array_cell(text, "material_id")
    gid = _array_cell(text, "global_element_id")
    if material is None:
        raise ValueError("NO CellData 'material_id'")
    if material.size != types.size:
        raise ValueError("material_id SIZE %d != CELL COUNT %d"
                         % (material.size, types.size))
    if gid is not None and gid.size != types.size:
        raise ValueError("global_element_id SIZE %d != CELL COUNT %d"
                         % (gid.size, types.size))
    tet = types == TET
    tets = conn[(offsets[tet] - 4)[:, None] + np.arange(4)[None, :]]
    material = material[tet]
    gid = gid[tet] if gid is not None else None
    return points, tets, temperature, material, gid


def read_result(path):
    path = latest(path)
    if path.endswith(".pvtu"):
        with open(path) as handle:
            sources = re.findall(r'<Piece\s+Source="([^"]+)"', handle.read())
        base = os.path.dirname(path)
        return [read_piece(os.path.join(base, s)) for s in sources], path
    return [read_piece(path)], path


def solid_energy(path, rho, cp, solid_id):
    pieces, path = read_result(path)
    cells = {}
    for ip, (points, tets, temperature, material, gid) in enumerate(pieces):
        keep = material == solid_id
        if not keep.any():
            continue
        tet = tets[keep]
        p = points[tet]
        vol = np.abs(np.einsum(
            "ij,ij->i",
            np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
            p[:, 3] - p[:, 0])) / 6.0
        tsum = temperature[tet].sum(axis=1)
        if gid is None:
            key = ["%d:%d" % (ip, i) for i in range(tet.shape[0])]
        else:
            key = gid[keep]
        for k, v, ts in zip(key, vol, tsum):
            cells[k] = (v, ts)       # GHOST COPIES OVERWRITE WITH THE SAME DATA
    if not cells:
        raise ValueError("NO SOLID TETRAHEDRA WITH material_id = %d" % solid_id)
    vol = np.array([v for v, _ in cells.values()])
    tsum = np.array([ts for _, ts in cells.values()])
    energy = rho * cp * np.sum(vol * tsum / 4.0)
    mean_t = np.sum(vol * tsum / 4.0) / np.sum(vol)
    return path, len(cells), np.sum(vol), mean_t, energy


def main():
    ap = argparse.ArgumentParser(description="CBS3D solid thermal energy")
    ap.add_argument("result", nargs="+")
    ap.add_argument("--rho", type=float, default=7850.0)
    ap.add_argument("--cp", type=float, default=500.0)
    ap.add_argument("--solid-id", type=int, default=1)
    ap.add_argument("--dt", type=float,
                    help="TIME BETWEEN THE FIRST AND LAST SNAPSHOT [s]")
    ap.add_argument("--qin", type=float,
                    help="HEAT INPUT [W], TO PRINT THE IMPLIED MEAN Q_out")
    args = ap.parse_args()

    name = "solid_id %d, rho %g kg/m3, cp %g J/kg/K" % (
        args.solid_id, args.rho, args.cp)
    print("SOLID   : %s" % name)
    print("%-58s %8s %14s %12s %16s"
          % ("FILE", "CELLS", "VOLUME [m3]", "MEAN T [K]", "E_s [J]"))
    rows = []
    for result in args.result:
        try:
            path, n, vol, mean_t, energy = solid_energy(
                result, args.rho, args.cp, args.solid_id)
        except ValueError as error:
            print("ERROR   : %s: %s" % (result, error))
            return 2
        rows.append((path, energy))
        print("%-58s %8d %14.6e %12.4f %16.8e"
              % (os.path.basename(path), n, vol, mean_t, energy))
    if len(rows) == 2:
        delta = rows[1][1] - rows[0][1]
        print()
        print("DELTA E : %.8e J" % delta)
        if args.dt is not None and args.dt > 0.0:
            storage = delta / args.dt
            print("STORAGE : %.6e W  (POSITIVE = SOLID WARMING)" % storage)
            if args.qin is not None:
                print("IMPLIED : <Q_out> = Q_in - storage = %.6e W"
                      % (args.qin - storage))
        else:
            print("RATE    : GIVE --dt TO CONVERT DELTA E TO W")
    return 0


if __name__ == "__main__":
    sys.exit(main())
