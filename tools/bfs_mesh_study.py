#!/usr/bin/env python3
# ================================
# FILE      : bfs_mesh_study.py
# PURPOSE   : GRID-INDEPENDENCE FIGURES FOR THE BACKWARD-FACING STEP: SEVERAL
#             MESHES (a, b, c, ... IN THE ORDER GIVEN, COARSE TO FINE)
#             AGAINST DENHAM ET AL. (1975).
#
#             <prefix>_profiles.png     u/U0 AT x/h = 0, 1, 2, 4, 6, 8 IN THE
#                                       LAYOUT OF LIU (2006) FIG. 4(b): DENHAM'S
#                                       DATA AND ONE CURVE PER MESH
#             <prefix>_convergence.png  VELOCITY RESIDUAL AGAINST TIME STEPS,
#                                       ONE CURVE PER MESH
#             TABLE (STDOUT)            ELEMENTS, ITERATIONS, x_r/h, RMS
#                                       AGAINST DENHAM PER STATION; FOR THREE
#                                       MESHES THE GRID CONVERGENCE INDEX OF
#                                       x_r/h (CELIK ET AL. 2008)
#
# USAGE     : python3 tools/bfs_mesh_study.py RUN_a RUN_b RUN_c [--out mesh_study]
#                 [--x-step 4] [--h 1] [--z 0.5] [--y-wall 0.02] [--target 1e-7]
#             RUN = THE run_<case>_<job> FOLDER OF EACH MESH
# ================================

from __future__ import print_function

import argparse
import glob
import os
import string
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from heat_transfer_coefficient import read_mesh  # noqa: E402
from probe_vtu import sample  # noqa: E402
from bfs_profiles import OFFSETS, read_ref  # noqa: E402
from plot_convergence import find_csv, read_column  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#c0392b", "#1baa7a", "#8e5bd0"]
STYLES = [(0, (6, 2, 1, 2)), (0, (5, 2)), "-", (0, (1, 1.5)), (0, (3, 1))]


def reattachment(X, U, tets, args):
    xs = args.x_step + args.h * np.linspace(0.01, 20.0, 2000)
    pts = np.c_[xs, np.full(len(xs), args.y_wall * args.h), np.full(len(xs), args.z)]
    uw = sample(X, U[:, 0:1], tets, pts)[:, 0]
    up = np.where((uw[:-1] < 0.0) & (uw[1:] >= 0.0))[0]
    if not len(up):
        return float("nan")
    i = up[-1]
    xr = xs[i] - uw[i] * (xs[i + 1] - xs[i]) / (uw[i + 1] - uw[i])
    return (xr - args.x_step) / args.h


def gci3(f, n):
    # CELIK ET AL. (2008): f, n = VALUES AND ELEMENT COUNTS, COARSE -> FINE
    from scipy.optimize import brentq
    f3, f2, f1 = f
    r21 = (n[2] / n[1]) ** (1.0 / 3.0)
    r32 = (n[1] / n[0]) ** (1.0 / 3.0)
    e32, e21 = f3 - f2, f2 - f1
    if e21 == 0.0 or e32 == 0.0:
        raise ValueError("NO CHANGE BETWEEN TWO MESHES")
    s = np.sign(e32 / e21)
    g = lambda p: p - abs(np.log(abs(e32 / e21)) + np.log((r21 ** p - s) / (r32 ** p - s))) / np.log(r21)
    p = brentq(g, 0.01, 20.0)
    fext = (r21 ** p * f1 - f2) / (r21 ** p - 1.0)
    return p, fext, 125.0 * abs((f1 - f2) / f1) / (r21 ** p - 1.0)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ref = os.path.join(here, "..", "cases", "backward_facing_step", "reference")
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", default="mesh_study")
    ap.add_argument("--x-step", type=float, default=4.0)
    ap.add_argument("--h", type=float, default=1.0)
    ap.add_argument("--z", type=float, default=0.5)
    ap.add_argument("--y-wall", type=float, default=0.02)
    ap.add_argument("--target", type=float, default=1.0e-7)
    ap.add_argument("--velocity", default="velocity")
    ap.add_argument("--ref", default=ref)
    args = ap.parse_args()

    ke, xh_e, ys_e, us_e = read_ref(os.path.join(args.ref, "liu2006_fig4b_exp.csv"))
    stations = sorted(set(zip(ke, xh_e)))

    meshes = []
    for i, run in enumerate(args.runs):
        files = sorted(glob.glob(os.path.join(run, "output", "*", "*_step_*.pvtu")))
        path, X, _T, U, _used, tets, _m, _bt, _bc = read_mesh(files[-1] if files else run, args.velocity, None)
        it, res = read_column(find_csv(run), "velocity_rel_max")
        letter = string.ascii_lowercase[i]
        meshes.append(dict(letter=letter, n=len(tets), X=X, U=U, tets=tets, it=it, res=res,
                           label="%s = %s elements" % (letter, format(len(tets), ","))))
        print("READ   : %s  %s" % (letter, path), file=sys.stderr)

    # ================================
    # TABLE
    # ================================
    print("%-4s %12s %11s %8s   rms vs Denham at x/h = %s" % (
        "mesh", "elements", "iterations", "x_r/h", ", ".join("%g" % s[1] for s in stations)))
    for m in meshes:
        m["xr"] = reattachment(m["X"], m["U"], m["tets"], args)
        rms = []
        for k, xh in stations:
            me = ke == k
            inner = (ys_e[me] > ys_e[me].min() + 1e-3) & (ys_e[me] < 3.0 - 1e-3)
            pts = np.c_[np.full(me.sum(), args.x_step + args.h * max(xh, 1e-6)), args.h * ys_e[me],
                        np.full(me.sum(), args.z)]
            u = sample(m["X"], m["U"][:, 0:1], m["tets"], pts)[:, 0]
            rms.append(np.sqrt(np.nanmean((u - us_e[me])[inner] ** 2)))
        print("%-4s %12s %11d %8.3f   %s" % (m["letter"], format(m["n"], ","), int(m["it"][-1]), m["xr"],
                                              "  ".join("%.4f" % r for r in rms)))
    if len(meshes) == 3:
        try:
            p, fext, g = gci3([m["xr"] for m in meshes], [float(m["n"]) for m in meshes])
            print("GCI    : x_r/h apparent order p = %.2f, extrapolated %.3f, GCI_fine = %.1f %%" % (p, fext, g))
        except Exception:
            print("GCI    : x_r/h not monotonically convergent over the three meshes")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "dejavuserif", "font.size": 10})

    # ================================
    # PROFILES: DENHAM + ONE CURVE PER MESH
    # ================================
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    for j, (k, xh) in enumerate(stations):
        off = OFFSETS[k]
        me = ke == k
        ax.plot(off + us_e[me], ys_e[me], "s", ms=3.2, mfc="none", mec="k", mew=0.7,
                label="Exp. (Denham et al.)" if j == 0 else None, zorder=5)
        y0 = args.h * (1.0 if xh <= 0.0 else 0.0)
        y_line = np.linspace(y0, 3.0 * args.h, 301)
        pts = np.c_[np.full(len(y_line), args.x_step + args.h * max(xh, 1e-6)), y_line, np.full(len(y_line), args.z)]
        for i, m in enumerate(meshes):
            u = sample(m["X"], m["U"][:, 0:1], m["tets"], pts)[:, 0]
            ax.plot(off + u, y_line / args.h, color=COLORS[i % len(COLORS)], ls=STYLES[i % len(STYLES)],
                    lw=1.2, label=m["label"] if j == 0 else None, zorder=3 + i)
        ax.text(off + 0.05, 3.06, "x/h = %g" % xh, fontsize=7.5, ha="left", va="bottom")
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 3)
    ax.set_xticks(range(0, 15, 2))
    ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
    ax.set_xlabel("Horizontal velocity")
    ax.set_ylabel("Vertical Distance")
    ax.set_title("Re = 3025", fontsize=10, pad=14)
    # LEGEND BELOW THE AXES, ONE ROW (CLEAR OF ALL PROFILES)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=len(meshes) + 1, fontsize=8,
              frameon=False, handlelength=2.2, columnspacing=1.6)
    fig.savefig(args.out + "_profiles.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ================================
    # CONVERGENCE HISTORIES
    # ================================
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    for i, m in enumerate(meshes):
        ax.semilogy(m["it"], m["res"], color=COLORS[i % len(COLORS)], ls=STYLES[i % len(STYLES)], lw=1.2,
                    label=m["label"])
    if args.target:
        ax.axhline(args.target, color="k", lw=0.7, ls=(0, (1, 2)),
                   label="convergence tolerance ($10^{%d}$)" % int(round(np.log10(args.target))))
    ax.set_xlabel("No. time steps")
    ax.set_ylabel(r"$\Vert e \Vert$ (velocity residual)")
    ax.set_title("Spalart-Allmaras model", fontsize=10)
    ax.set_xlim(left=0)
    ax.grid(True, which="major", color="#d8d8d4", lw=0.5)
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    fig.savefig(args.out + "_convergence.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("WROTE  : %s_profiles.png, %s_convergence.png" % (args.out, args.out))


if __name__ == "__main__":
    main()
