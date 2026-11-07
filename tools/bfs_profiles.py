#!/usr/bin/env python3
# ================================
# FILE      : bfs_profiles.py
# PURPOSE   : BACKWARD-FACING STEP, Re = 3025: COMPARE A CBS3D RESULT WITH
#             LIU'S SA SOLUTION AND DENHAM'S MEASUREMENTS
#             (cases/backward_facing_step/reference, FROM
#             tools/digitise_liu2006_bfs.py).
#
#             AT EACH STATION x/h = 0, 1, 2, 4, 6, 8 DOWNSTREAM OF THE STEP,
#             u/U0 IS SAMPLED (P1 INTERPOLATION, MID-SPAN z) AT THE y OF EVERY
#             REFERENCE POINT, AND
#               rms_liu, max_liu   AGAINST LIU'S SA CURVE
#               rms_exp            AGAINST DENHAM'S DATA (WALL POINTS EXCLUDED)
#               Q                  INT u dy OVER THE SECTION (LIU SA: 1.93-1.96)
#             ARE PRINTED. THE REATTACHMENT LENGTH x_r/h IS THE LAST SIGN
#             CHANGE (- TO +) OF u ALONG THE FLOOR AT HEIGHT --y-wall, WITHIN
#             20 h OF THE STEP (DENHAM: ABOUT 6 h).
#
# USAGE     : python3 tools/bfs_profiles.py RESULT [--x-step 4] [--h 1]
#                 [--z 0.5] [--y-wall 0.02] [--velocity velocity]
#                 [--csv out.csv] [--plot out.png]
#             RESULT = .pvtu / .vtu / FOLDER (LATEST .pvtu UNDER IT)
#             --plot WRITES A FIGURE IN THE LAYOUT OF LIU (2006) FIG. 4(b)
#
#             WALL LAW (--wall-law FILE.png): u+ = u/u_tau AGAINST
#             y+ = y u_tau/nu AT ATTACHED STATIONS (--wall-stations, DEFAULT
#             x/h = -2 FLOOR AND ROOF, 12, 20, 30 FLOOR, 30 ROOF; NOT INSIDE
#             THE RECIRCULATION, WHERE u_tau -> 0). u_tau FROM THE EXACT P1
#             WALL GRADIENT: nu_t = 0 AT THE WALL, SO tau_w = mu du/dn, AND u
#             IS LINEAR IN THE FIRST ELEMENT: du/dn = u(delta)/delta. ALSO
#             PRINTS Re_tau (HALF CHANNEL) AND y+ OF THE NEAREST OFF-WALL NODE.
#             REFERENCE LINES AS LIU: u+ = y+, u+ = (1/0.418) ln y+ + 5.5.
# ================================

from __future__ import print_function

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from heat_transfer_coefficient import read_mesh  # noqa: E402
from probe_vtu import sample  # noqa: E402

OFFSETS = [0.0, 2.0, 4.0, 6.0, 8.0, 12.0]       # PLOT SHIFTS OF LIU FIG. 4(b)


def read_ref(path):
    # "#" COMMENT LINES, THEN THE HEADER station,x_over_h,y,u, THEN DATA
    rows = [ln for ln in open(path) if ln.strip() and not ln.startswith("#")]
    d = np.array([[float(v) for v in ln.split(",")] for ln in rows[1:]])
    return d[:, 0].astype(int), d[:, 1], d[:, 2], d[:, 3]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ref = os.path.join(here, "..", "cases", "backward_facing_step", "reference")
    ap = argparse.ArgumentParser()
    ap.add_argument("result")
    ap.add_argument("--x-step", type=float, default=4.0, help="x OF THE STEP FACE IN THE MESH")
    ap.add_argument("--h", type=float, default=1.0, help="STEP HEIGHT")
    ap.add_argument("--z", type=float, default=0.5, help="SAMPLING PLANE (MID-SPAN)")
    ap.add_argument("--y-wall", type=float, default=0.02, help="HEIGHT OF THE REATTACHMENT LINE [h]")
    ap.add_argument("--velocity", default="velocity")
    ap.add_argument("--ref", default=ref)
    ap.add_argument("--csv", default=None)
    ap.add_argument("--plot", default=None)
    ap.add_argument("--wall-law", default=None, help="FIGURE FILE FOR u+ AGAINST y+")
    ap.add_argument("--wall-stations", default="-2:floor,-2:roof,12:floor,20:floor,30:floor,30:roof")
    ap.add_argument("--mu", type=float, default=1.0 / 3025.0, help="MOLECULAR VISCOSITY (rho = 1)")
    args = ap.parse_args()

    path, X, _T, U, _used, tets, _mats, _bt, _bc = read_mesh(args.result, args.velocity, None)
    print("READ   : %s (%d nodes, %d tetrahedra)" % (path, len(X) , len(tets)), file=sys.stderr)

    ks, xh_s, ys_s, us_s = read_ref(os.path.join(args.ref, "liu2006_fig4b_sa.csv"))
    ke, xh_e, ys_e, us_e = read_ref(os.path.join(args.ref, "liu2006_fig4b_exp.csv"))

    def u_at(xh, y):
        # A STATION ON THE STEP FACE (x/h = 0) IS SAMPLED A HAIR DOWNSTREAM
        x = args.x_step + args.h * max(xh, 1.0e-6)
        pts = np.c_[np.full(len(y), x), args.h * np.asarray(y), np.full(len(y), args.z)]
        return sample(X, U[:, 0:1], tets, pts)[:, 0]

    rows = []
    print("%-6s %8s %8s %8s %8s %8s %6s" % ("x/h", "rms_liu", "max_liu", "rms_exp", "Q_cbs", "Q_liu", "n_nan"))
    for k in range(len(OFFSETS)):
        ms, me = ks == k, ke == k
        xh = xh_s[ms][0]
        u_s = u_at(xh, ys_s[ms])
        u_e = u_at(xh, ys_e[me])
        inner = (ys_e[me] > ys_e[me].min() + 1e-3) & (ys_e[me] < 3.0 - 1e-3)
        dl = u_s - us_s[ms]
        de = (u_e - us_e[me])[inner]
        y_line = np.linspace(ys_s[ms].min(), ys_s[ms].max(), 301)
        q_cbs = np.trapezoid(np.nan_to_num(u_at(xh, y_line)), y_line)
        q_liu = np.trapezoid(us_s[ms], ys_s[ms])
        nn = int(np.isnan(u_s).sum() + np.isnan(u_e).sum())
        print("%-6g %8.4f %8.4f %8.4f %8.4f %8.4f %6d" % (
            xh, np.sqrt(np.nanmean(dl ** 2)), np.nanmax(np.abs(dl)), np.sqrt(np.nanmean(de ** 2)),
            q_cbs, q_liu, nn))
        for y, a, b in zip(ys_s[ms], u_s, us_s[ms]):
            rows.append((xh, "liu_sa", y, a, b))
        for y, a, b in zip(ys_e[me], u_e, us_e[me]):
            rows.append((xh, "denham", y, a, b))

    # ================================
    # REATTACHMENT: LAST - TO + SIGN CHANGE OF u ALONG THE FLOOR
    # ================================
    xs = args.x_step + args.h * np.linspace(0.01, 20.0, 2000)
    pts = np.c_[xs, np.full(len(xs), args.y_wall * args.h), np.full(len(xs), args.z)]
    uw = sample(X, U[:, 0:1], tets, pts)[:, 0]
    up = np.where((uw[:-1] < 0.0) & (uw[1:] >= 0.0))[0]
    if len(up):
        i = up[-1]
        xr = xs[i] - uw[i] * (xs[i + 1] - xs[i]) / (uw[i + 1] - uw[i])
        print("x_r/h  : %.3f   (u = 0 AT y = %g h; DENHAM ~6)" % ((xr - args.x_step) / args.h, args.y_wall))
    else:
        print("x_r/h  : NO REVERSED FLOW FOUND AT y = %g h" % args.y_wall)

    if args.wall_law:
        wall_law(args, X, U, tets)

    if args.csv:
        with open(args.csv, "w") as f:
            f.write("x_over_h,reference,y,u_cbs,u_reference\n")
            for r in rows:
                f.write("%g,%s,%.5f,%.6f,%.6f\n" % r)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # ================================
        # LAYOUT OF LIU (2006) FIG. 4(b) / THESIS FIG. 8.3: EACH PROFILE
        # SHIFTED BY ITS OFFSET, AXES 0..14 x 0..3, TITLE Re = 3025
        # ================================
        plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "dejavuserif", "font.size": 10})
        fig, ax = plt.subplots(figsize=(7.2, 3.2))
        for k, off in enumerate(OFFSETS):
            ms, me = ks == k, ke == k
            y_line = np.linspace(ys_s[ms].min(), 3.0, 301)
            ax.plot(off + us_e[me], ys_e[me], "s", ms=3.2, mfc="none", mec="k", mew=0.7,
                    label="Exp. (Denham et al.)" if k == 0 else None, zorder=3)
            ax.plot(off + us_s[ms], ys_s[ms], color="k", lw=0.9, ls=(0, (4, 2)),
                    label="SA, Liu (2006)" if k == 0 else None, zorder=2)
            ax.plot(off + u_at(xh_s[ms][0], y_line), y_line, color="#c0392b", lw=1.3,
                    label="SA, CBS3D" if k == 0 else None, zorder=4)
            ax.text(off + 0.05, 3.06, "x/h = %g" % xh_s[ms][0], fontsize=7.5, ha="left", va="bottom")
        ax.set_xlim(0, 14)
        ax.set_ylim(0, 3)
        ax.set_xticks(range(0, 15, 2))
        ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
        ax.set_xlabel("Horizontal velocity")
        ax.set_ylabel("Vertical Distance")
        ax.set_title("Re = 3025", fontsize=10, pad=14)
        # LEGEND IN THE EMPTY GAP BETWEEN THE x/h = 6 AND 8 PROFILES
        ax.legend(loc="center", bbox_to_anchor=(10.95 / 14.0, 0.5), fontsize=7.5, frameon=False,
                  handlelength=1.8, borderaxespad=0.0)
        fig.savefig(args.plot, dpi=300, bbox_inches="tight")


def locate(X, tets, p):
    # INDEX OF THE TETRAHEDRON CONTAINING p (BARYCENTRIC TEST), OR None
    P = X[tets]
    cand = np.where(np.all(P.min(axis=1) <= p + 1e-12, axis=1) & np.all(P.max(axis=1) >= p - 1e-12, axis=1))[0]
    for e in cand:
        A = np.vstack([np.ones(4), X[tets[e]].T])
        w = np.linalg.solve(A, np.r_[1.0, p])
        if w.min() >= -1e-9:
            return e
    return None


def wall_law(args, X, U, tets):
    # ================================
    # u+ AGAINST y+ AT ATTACHED STATIONS
    # ================================
    nu = args.mu
    curves = []
    print("\n%-6s %-6s %10s %8s %8s %10s" % ("x/h", "wall", "u_tau", "Re_tau", "y1+_elem", "note"))
    for item in args.wall_stations.split(","):
        xs, wall = item.split(":")
        xh = float(xs)
        x = args.x_step + args.h * xh
        floor = args.h if x < args.x_step else 0.0
        roof = 3.0 * args.h
        y_w, sgn = (floor, 1.0) if wall == "floor" else (roof, -1.0)
        half = 0.5 * (roof - floor)

        def u_n(d):
            pts = np.c_[np.full(len(d), x), y_w + sgn * np.asarray(d), np.full(len(d), args.z)]
            return sample(X, U[:, 0:1], tets, pts)[:, 0]

        delta = 1.0e-5 * args.h
        dudn = u_n(np.array([delta]))[0] / delta
        # FIRST-ELEMENT HEIGHT: THE WALL-ADJACENT TETRAHEDRON CONTAINING THE
        # POINT delta ABOVE THE WALL; ITS HIGHEST NODE ABOVE THE WALL
        e = locate(X, tets, np.array([x, y_w + sgn * delta, args.z]))
        dist = sgn * (X[tets[e], 1] - y_w) if e is not None else np.array([])
        if not np.isfinite(dudn) or dudn <= 0.0:
            print("%-6g %-6s %10s %8s %8s %10s" % (xh, wall, "-", "-", "-", "REVERSED"))
            continue
        u_tau = np.sqrt(nu * dudn)
        y1p = dist.max() * u_tau / nu if len(dist) else float("nan")
        yp = np.logspace(-0.3, np.log10(half * u_tau / nu), 240)
        up = u_n(yp * nu / u_tau) / u_tau
        curves.append((xh, wall, yp, up))
        print("%-6g %-6s %10.5f %8.1f %8.2f %10s" % (xh, wall, u_tau, half * u_tau / nu, y1p, ""))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "dejavuserif", "font.size": 10})
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    yl = np.logspace(0, 3, 200)
    ax.plot(yl[yl < 15], yl[yl < 15], color="k", lw=0.8, ls=(0, (1, 2)), label="$u^+ = y^+$")
    ax.plot(yl[yl > 5], np.log(yl[yl > 5]) / 0.418 + 5.5, color="k", lw=0.8, ls=(0, (5, 3)),
            label="$u^+ = (1/0.418)\\ln y^+ + 5.5$")
    colors = ["#2a78d6", "#eb6834", "#1baa7a", "#8e5bd0", "#c0392b", "#7a7a72"]
    for i, (xh, wall, yp, up) in enumerate(curves):
        ax.plot(yp, up, color=colors[i % len(colors)], lw=1.4, label="x/h = %g, %s" % (xh, wall))
    ax.set_xscale("log")
    ax.set_xlim(1, 1000)
    ax.set_ylim(0, 25)
    ax.set_xlabel("$y^+$")
    ax.set_ylabel("$u^+$")
    ax.set_title("Re = 3025; Spalart-Allmaras model", fontsize=10)
    ax.grid(True, which="major", color="#d8d8d4", lw=0.5)
    ax.legend(loc="upper left", fontsize=7.5, frameon=False)
    fig.savefig(args.wall_law, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
