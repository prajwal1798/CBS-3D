#!/usr/bin/env python3
# ================================
# FILE      : check_nuclear_slab.py
# PURPOSE   : COMPARE THE CBS3D STEADY TEMPERATURE OF THE nuclear_heating_slab
#             CASE WITH THE EXACT 1D SOLUTION
#               k T'' = -A exp(-(L-y)/lam),  T(0) = 0,  k T'(L) = q
#               T(y)  = (A/k)[lam y - lam^2 (exp(-(L-y)/lam) - exp(-L/lam))]
#                       + q y / k
#             PASS WHEN THE LARGEST NODAL ERROR IS BELOW --tol OF T(L).
#             --compare OTHER ALSO REQUIRES THE SAME NODAL TEMPERATURES AS
#             ANOTHER RUN (E.G. SERIAL VS MPI) TO WITHIN 1e-6 OF T(L).
#
# USAGE     : python3 tools/check_nuclear_slab.py RESULT [--A 1.6e6]
#                 [--lam 0.05] [--q 1e5] [--k 28.3] [--L 0.10] [--W 0.01]
#                 [--tol 0.01] [--compare OTHER_RESULT]
#             RESULT = A .vtu, A .pvtu OR A FOLDER (LATEST FILE IS USED)
#
# STATUS    : NUMPY ONLY. CBS3D ASCII VTU.
# ================================

from __future__ import print_function

import argparse
import os
import re
import sys

import numpy as np


# ================================
# CBS3D ASCII VTU READER
# ================================
def _array(text, name):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                  text, re.S)
    if m is None:
        raise ValueError("NO DataArray '%s'" % name)
    return np.array(m.group(1).split(), dtype=float)


def read_piece(path):
    """POINTS (n,3) AND NODAL TEMPERATURE (n,) OF ONE PIECE."""
    with open(path) as handle:
        text = handle.read()
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    points = np.array(m.group(1).split(), dtype=float).reshape(-1, 3)
    return points, _array(text, "temperature")


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


def read_result(path):
    path = latest(path)
    print("READ   :", path)
    if path.endswith(".pvtu"):
        with open(path) as handle:
            sources = re.findall(r'<Piece\s+Source="([^"]+)"', handle.read())
        base = os.path.dirname(path)
        pieces = [read_piece(os.path.join(base, s)) for s in sources]
        return (np.vstack([p[0] for p in pieces]),
                np.concatenate([p[1] for p in pieces]))
    return read_piece(path)


# ================================
# EXACT SOLUTION
# ================================
def exact(y, a, lam, q, k, length):
    return (a / k) * (lam * y - lam ** 2 * (np.exp(-(length - y) / lam)
                                            - np.exp(-length / lam))) + q * y / k


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("result")
    ap.add_argument("--A", type=float, default=1.6e6)
    ap.add_argument("--lam", type=float, default=0.05)
    ap.add_argument("--q", type=float, default=1.0e5)
    ap.add_argument("--k", type=float, default=28.3)
    ap.add_argument("--L", type=float, default=0.10)
    ap.add_argument("--W", type=float, default=0.01)
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--compare", help="ANOTHER RESULT THAT MUST AGREE NODE BY NODE")
    args = ap.parse_args()

    points, temp = read_result(args.result)
    bar = points[:, 0] <= args.W * (1.0 + 1e-6)          # SOLID BAR ONLY
    y = points[bar, 1]
    t_num = temp[bar]
    t_ref = exact(y, args.A, args.lam, args.q, args.k, args.L)
    t_top = exact(np.array([args.L]), args.A, args.lam, args.q, args.k, args.L)[0]

    err = np.abs(t_num - t_ref)
    worst = int(np.argmax(err))
    rel = err[worst] / abs(t_top)

    print("NODES  : %d in the solid bar" % y.size)
    print("T(L)   : exact %.4f K   CBS3D max %.4f K" % (t_top, t_num.max()))
    print("ERROR  : max |T - T_exact| = %.4e K at y = %.4f m  (%.3e of T(L))"
          % (err[worst], y[worst], rel))
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        near = np.argmin(np.abs(y - frac * args.L))
        print("   y = %.3f  T = %10.4f  exact = %10.4f"
              % (y[near], t_num[near], t_ref[near]))
    ok = rel <= args.tol
    if args.compare:
        other_pts, other_t = read_result(args.compare)
        lookup = {}
        for p, value in zip(other_pts, other_t):
            lookup[tuple(np.round(p, 9))] = value
        worst_diff = 0.0
        for p, value in zip(points, temp):
            key = tuple(np.round(p, 9))
            if key not in lookup:
                print("COMPARE: NODE %s MISSING IN %s" % (key, args.compare))
                return 1
            worst_diff = max(worst_diff, abs(value - lookup[key]))
        same = worst_diff <= 1.0e-6 * abs(t_top)
        print("COMPARE: max |T - T_other| = %.3e K (%s)"
              % (worst_diff, "SAME" if same else "DIFFERENT"))
        ok = ok and same
    if ok:
        print("RESULT : PASS (tolerance %.1e)" % args.tol)
        return 0
    print("RESULT : FAIL (tolerance %.1e)" % args.tol)
    return 1


if __name__ == "__main__":
    sys.exit(main())
