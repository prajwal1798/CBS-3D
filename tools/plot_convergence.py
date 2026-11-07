#!/usr/bin/env python3
# ================================
# FILE      : plot_convergence.py
# PURPOSE   : STEADY-STATE CONVERGENCE HISTORIES OF ONE OR MORE CBS3D RUNS IN
#             THE LAYOUT OF LIU (2005) FIG. 8.4 / NITHIARASU & LIU (2006)
#             FIG. 3: VELOCITY RESIDUAL (LOG SCALE) AGAINST TIME STEPS, ONE
#             CURVE PER RUN.
#
#             THE CURVE IS velocity_rel_max OF THE RESIDUAL CSV: THE LARGEST
#             OF THE RELATIVE CHANGES OF u, v, w IN ONE ITERATION, I.E. THE
#             QUANTITY THE STOP CRITERION TESTS. LIU PLOTS THE L2 NORM OF HIS
#             VELOCITY RESIDUAL; THE NORMALISATIONS DIFFER, SO COMPARE SHAPES
#             AND ITERATION COUNTS, NOT LEVELS. --column CHOOSES ANOTHER COLUMN
#             (E.G. u_rel, sa_rel).
#
# USAGE     : python3 tools/plot_convergence.py RUN_DIR[=LABEL] [RUN_DIR[=LABEL] ...]
#                 [--out convergence.png] [--column velocity_rel_max]
#                 [--title "Spalart-Allmaras model"] [--target 1e-7]
#             RUN_DIR = run_<case>_<job> FOLDER, output/<case> FOLDER, OR THE
#             *_residuals.csv FILE ITSELF
# ================================

from __future__ import print_function

import argparse
import glob
import os
import sys

import numpy as np


def find_csv(path):
    if os.path.isfile(path):
        return path
    for pattern in ("output/*/*_residuals.csv", "*_residuals.csv"):
        hits = sorted(glob.glob(os.path.join(path, pattern)))
        if hits:
            return hits[0]
    sys.exit("NO *_residuals.csv UNDER %s" % path)


def read_column(path, name):
    with open(path) as f:
        head = f.readline().strip().split(",")
    if name not in head:
        sys.exit("NO COLUMN %s IN %s (HAS: %s)" % (name, path, ", ".join(head)))
    d = np.genfromtxt(path, delimiter=",", skip_header=1, usecols=(head.index("iteration"), head.index(name)))
    d = d[np.isfinite(d).all(axis=1) & (d[:, 1] > 0)]
    return d[:, 0], d[:, 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", default="convergence.png")
    ap.add_argument("--column", default="velocity_rel_max")
    ap.add_argument("--title", default="Spalart-Allmaras model")
    ap.add_argument("--target", type=float, default=None, help="DRAW THE STOP TOLERANCE")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "dejavuserif", "font.size": 10})
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    colors = ["#c0392b", "#2a78d6", "#1baa7a", "#8e5bd0", "#eb6834", "#55554f"]
    styles = ["-", (0, (5, 2)), (0, (1, 1.5)), (0, (6, 2, 1, 2)), "-", (0, (3, 1))]

    for i, item in enumerate(args.runs):
        path, label = (item.split("=", 1) + [None])[:2]
        csv = find_csv(path)
        it, r = read_column(csv, args.column)
        label = label or os.path.basename(os.path.normpath(path))
        ax.semilogy(it, r, color=colors[i % len(colors)], ls=styles[i % len(styles)], lw=1.2, label=label)
        print("%-40s %8d iterations, final %s = %.3e" % (label, int(it[-1]), args.column, r[-1]))

    if args.target:
        ax.axhline(args.target, color="k", lw=0.6, ls=(0, (1, 2)))
    ax.set_xlabel("No. time steps")
    ax.set_ylabel(r"$\Vert e \Vert$ (velocity residual)" if args.column == "velocity_rel_max" else args.column)
    ax.set_title(args.title, fontsize=10)
    ax.set_xlim(left=0)
    ax.grid(True, which="major", color="#d8d8d4", lw=0.5)
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    print("WROTE  : %s" % args.out)


if __name__ == "__main__":
    main()
