#!/usr/bin/env python3
# ================================
# FILE      : probe_vtu.py
# PURPOSE   : TEMPERATURE (AND SPEED) AT THE PROBE POINTS OF A PROBE FILE,
#             READ FROM SAVED CBS3D RESULTS, ONE COLUMN PER RESULT. COMPARES
#             RUNS THAT HAD NO IN-SOLVER PROBES (E.G. STEADY ENERGY RUNS)
#             WITH THE PROBE TRACES OF A TRANSIENT RUN.
#
#             EACH POINT IS LOCATED IN ITS TETRAHEDRON AND THE NODAL VALUES
#             ARE INTERPOLATED WITH THE BARYCENTRIC (P1) WEIGHTS, THE SAME
#             RULE AS THE SOLVER'S PointSampler.
#
# USAGE     : python3 tools/probe_vtu.py RESULT [RESULT ...]
#                 [--probes cases/blanket_slice/blanket_probes.txt]
#                 [--velocity velocity] [--csv out.csv]
#             RESULT = .pvtu / .vtu / FOLDER (LATEST FILE UNDER IT)
#             RESULT:mean  READS mean_temperature AND mean_velocity (THE TIME
#                          AVERAGES OF A RUN WITH CBS3D_AVERAGE_START_TIME)
#                          INSTEAD OF THE INSTANTANEOUS FIELDS
#
# STATUS    : NUMPY ONLY (READER SHARED WITH heat_transfer_coefficient.py).
# ================================

from __future__ import print_function

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from heat_transfer_coefficient import read_mesh  # noqa: E402


def read_probes(path):
    names, xyz = [], []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        p = s.split()
        names.append(p[0])
        xyz.append([float(v) for v in p[1:4]])
    return names, np.array(xyz)


def sample(X, field, tets, points):
    lo = X[tets].min(axis=1)
    hi = X[tets].max(axis=1)
    size = (hi - lo).max(axis=1)
    out = np.full((len(points),) + field.shape[1:], np.nan)
    for i, p in enumerate(points):
        tol = 1e-9 * size
        cand = np.where(np.all(lo <= p + tol[:, None], axis=1) & np.all(hi >= p - tol[:, None], axis=1))[0]
        for e in cand:
            n = tets[e]
            A = np.vstack([np.ones(4), X[n].T])
            w = np.linalg.solve(A, np.r_[1.0, p])
            if w.min() >= -1e-9:
                out[i] = np.tensordot(w, field[n], axes=1)
                break
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="+")
    ap.add_argument("--probes", default=os.path.join(here, "..", "cases", "blanket_slice", "blanket_probes.txt"))
    ap.add_argument("--velocity", default="velocity")
    ap.add_argument("--csv", default=None)
    args = ap.parse_args()

    names, pts = read_probes(args.probes)
    cols = []
    label = []
    for r in args.results:
        tname, vname, tag = "temperature", args.velocity, ""
        if r.endswith(":mean"):
            r, tname, vname, tag = r[:-5], "mean_temperature", "mean_velocity", "_mean"
        path, X, T, U, _used, tets, _mats, _bt, _bc = read_mesh(r, vname, tname)
        both = np.c_[T, np.linalg.norm(U, axis=1)]
        cols.append((path, sample(X, both, tets, pts)))
        label.append((os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(path)))) or path) + tag)
        print("READ  : %s (%s, %s)" % (path, tname, vname), file=sys.stderr)

    print("%-16s" % "probe" + "".join("%22s" % ("T [K] " + l[-14:]) for l in label)
          + "".join("%20s" % ("|u| " + l[-14:]) for l in label))
    for i, n in enumerate(names):
        print("%-16s" % n + "".join("%22.3f" % c[1][i, 0] for c in cols)
              + "".join("%20.5f" % c[1][i, 1] for c in cols))
    if args.csv:
        with open(args.csv, "w") as f:
            f.write("probe," + ",".join("T_%s,speed_%s" % (l, l) for l in label) + "\n")
            for i, n in enumerate(names):
                f.write(n + "," + ",".join("%.6f,%.6f" % (c[1][i, 0], c[1][i, 1]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
