#!/usr/bin/env python3
# ================================
# FILE      : plot_ldc_centerlines.py
# PURPOSE   : LID-DRIVEN CAVITY CENTRELINE PROFILES FROM A CBS3D RESULT,
#             AGAINST GHIA, GHIA AND SHIN (1982), IN A PUBLICATION STYLE:
#               u ALONG THE VERTICAL CENTRELINE   (x = xc, z = zc)
#               v ALONG THE HORIZONTAL CENTRELINE (y = yc, z = zc)
#
# USAGE     : python plot_ldc_centerlines.py RESULT [OPTIONS]
#
#   RESULT    A .pvtu (MPI RUN), A .vtu (SERIAL RUN) OR AN output/<case>
#             FOLDER (THE LAST STEP IN IT IS USED).
#
# OPTIONS   : --previous RESULT   AN EARLIER STEP OF THE SAME RUN. PRINTS HOW
#                                 MUCH THE PROFILES STILL CHANGE (CONVERGENCE).
#             --re 100            GHIA DATA SET: 100 OR 400 (MUST MATCH THE .par!)
#             --ref-u CSV         OTHER REFERENCE u PROFILE: y,u PER LINE
#             --ref-v CSV         OTHER REFERENCE v PROFILE: x,v PER LINE
#             --ref-label TEXT    LEGEND TEXT FOR --ref-u / --ref-v
#             --label TEXT        LEGEND TEXT FOR THE SOLUTION (DEFAULT CBS3D)
#             --z Z               DEPTH OF THE PROFILES (DEFAULT: MID-DEPTH)
#             --out PREFIX        OUTPUT NAME (DEFAULT ldc_re<Re>)
#             --samples N         POINTS PER PROFILE (DEFAULT 401)
#
# OUTPUT    : PREFIX_u.png/.pdf, PREFIX_v.png/.pdf  (600 dpi PNG, VECTOR PDF)
#             PREFIX_profiles.csv                   (THE SAMPLED PROFILES)
#
# NEEDS     : numpy, matplotlib. NO VTK / PARAVIEW: READS CBS3D ASCII VTU
#             DIRECTLY AND INTERPOLATES INSIDE THE TETRAHEDRA (EXACT P1 VALUES).
#
# NOTE      : GHIA ET AL. IS A 2D SOLUTION. A 3D CUBE (SIDE WALLS) GIVES A
#             WEAKER VORTEX ON THE MID-PLANE, SO SMALLER |v| PEAKS ARE EXPECTED.
# ================================

from __future__ import print_function

import argparse
import os
import re
import sys

import numpy as np

# ================================
# GHIA, GHIA AND SHIN (1982), J. COMPUT. PHYS. 48, 387-411, TABLES I AND II
# ================================
GHIA = {
    100: {
        "u": np.array([
            [1.0000, 1.00000], [0.9766, 0.84123], [0.9688, 0.78871],
            [0.9609, 0.73722], [0.9531, 0.68717], [0.8516, 0.23151],
            [0.7344, 0.00332], [0.6172, -0.13641], [0.5000, -0.20581],
            [0.4531, -0.21090], [0.2813, -0.15662], [0.1719, -0.10150],
            [0.1016, -0.06434], [0.0703, -0.04775], [0.0625, -0.04192],
            [0.0547, -0.03717], [0.0000, 0.00000]]),
        "v": np.array([
            [1.0000, 0.00000], [0.9688, -0.05906], [0.9609, -0.07391],
            [0.9531, -0.08864], [0.9453, -0.10313], [0.9063, -0.16914],
            [0.8594, -0.22445], [0.8047, -0.24533], [0.5000, 0.05454],
            [0.2344, 0.17527], [0.2266, 0.17507], [0.1563, 0.16077],
            [0.0938, 0.12317], [0.0781, 0.10890], [0.0703, 0.10091],
            [0.0625, 0.09233], [0.0000, 0.00000]]),
    },
    400: {
        "u": np.array([
            [1.0000, 1.00000], [0.9766, 0.75837], [0.9688, 0.68439],
            [0.9609, 0.61756], [0.9531, 0.55892], [0.8516, 0.29093],
            [0.7344, 0.16256], [0.6172, 0.02135], [0.5000, -0.11477],
            [0.4531, -0.17119], [0.2813, -0.32726], [0.1719, -0.24299],
            [0.1016, -0.14612], [0.0703, -0.10338], [0.0625, -0.09266],
            [0.0547, -0.08186], [0.0000, 0.00000]]),
        # x = 0.9063 (-0.23827) IS WIDELY REGARDED AS A MISPRINT IN THE PAPER;
        # IT IS KEPT AS PUBLISHED.
        "v": np.array([
            [1.0000, 0.00000], [0.9688, -0.12146], [0.9609, -0.15663],
            [0.9531, -0.19254], [0.9453, -0.22847], [0.9063, -0.23827],
            [0.8594, -0.44993], [0.8047, -0.38598], [0.5000, 0.05186],
            [0.2344, 0.30174], [0.2266, 0.30203], [0.1563, 0.28124],
            [0.0938, 0.22965], [0.0781, 0.20920], [0.0703, 0.19713],
            [0.0625, 0.18360], [0.0000, 0.00000]]),
    },
}

# ================================
# STYLE
# ================================
BACKGROUND = "#F3F0FA"      # VERY MILD VIOLET
LINE_COLOUR = "#1F4E9C"     # SOLUTION
MARK_COLOUR = "#E0661B"     # REFERENCE
FONTS = ["Arial", "Liberation Sans", "Helvetica", "Arimo", "DejaVu Sans"]


def set_style(plt):
    from matplotlib import font_manager
    # WSL / WINDOWS: USE THE REAL ARIAL FROM THE WINDOWS FONT FOLDER.
    for folder in ("/mnt/c/Windows/Fonts", "C:/Windows/Fonts"):
        for name in ("arial.ttf", "ariali.ttf", "arialbd.ttf", "arialbi.ttf"):
            path = os.path.join(folder, name)
            if os.path.isfile(path):
                try:
                    font_manager.fontManager.addfont(path)
                except Exception:
                    pass
    have = {f.name for f in font_manager.fontManager.ttflist}
    family = next((f for f in FONTS if f in have), "DejaVu Sans")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [family],
        # LaTeX-STYLE MATHS ($u$, $y$) SET IN THE SAME SANS FONT.
        "mathtext.fontset": "custom",
        "mathtext.rm": family,
        "mathtext.it": family + ":italic",
        "mathtext.bf": family + ":bold",
        "mathtext.sf": family,
        "font.size": 9,
        "axes.labelsize": 10,
        "legend.fontsize": 8.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "black",
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.size": 2.0,
        "ytick.minor.size": 2.0,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
    })
    return family


def style_axes(ax):
    from matplotlib.ticker import AutoMinorLocator
    ax.set_facecolor(BACKGROUND)
    ax.set_axisbelow(True)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(which="major", color="white", linewidth=0.9)
    ax.grid(which="minor", color="white", linewidth=0.45)
    ax.tick_params(which="both", top=True, right=True)
    for side in ax.spines.values():
        side.set_color("black")


def legend_above(ax):
    # OUTSIDE THE AXES, ABOVE THE FRAME: IT CAN NEVER HIDE A CURVE OR MARKER.
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2,
              frameon=False, handlelength=2.2, columnspacing=1.6,
              borderaxespad=0.2)


# ================================
# CBS3D ASCII VTU READER
# ================================
def _array(text, name):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                  text, re.S)
    if m is None:
        raise ValueError("NO DataArray '%s'" % name)
    return np.array(m.group(1).split(), dtype=float)


def read_vtu(path):
    """POINTS (n,3), TETRAHEDRA (m,4) AND NODAL VELOCITY (n,3) OF ONE PIECE."""
    with open(path) as handle:
        text = handle.read()
    if 'format="ascii"' not in text:
        raise ValueError("%s: ONLY ASCII VTU IS SUPPORTED" % path)
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    points = np.array(m.group(1).split(), dtype=float).reshape(-1, 3)
    conn = _array(text, "connectivity").astype(np.int64)
    offsets = _array(text, "offsets").astype(np.int64)
    types = _array(text, "types").astype(np.int64)
    tet = types == 10                        # VTK_TETRA; SKIPS BOUNDARY FACES
    start = offsets[tet] - 4
    tets = conn[start[:, None] + np.arange(4)[None, :]]
    velocity = _array(text, "velocity").reshape(-1, 3)
    return points, tets, velocity


def read_result(path):
    """LIST OF (POINTS, TETS, VELOCITY) PIECES."""
    if os.path.isdir(path):
        files = sorted(f for f in os.listdir(path) if f.endswith(".pvtu"))
        if not files:
            files = sorted(f for f in os.listdir(path)
                           if f.endswith(".vtu") and "_rank_" not in f
                           and "boundary" not in f)
        if not files:
            raise ValueError("NO .pvtu OR .vtu IN %s" % path)
        path = os.path.join(path, files[-1])
    print("READ   :", path)
    if path.endswith(".pvtu"):
        with open(path) as handle:
            sources = re.findall(r'<Piece\s+Source="([^"]+)"', handle.read())
        base = os.path.dirname(path)
        return [read_vtu(os.path.join(base, s)) for s in sources], path
    return [read_vtu(path)], path


# ================================
# P1 INTERPOLATION ALONG A STRAIGHT LINE
# ================================
def sample_line(pieces, fixed, axis, coords):
    """VELOCITY AT POINTS WHERE coordinate[axis] = coords AND THE OTHER TWO
    COORDINATES ARE fixed = {axis: value}. NaN WHERE NO TETRAHEDRON IS FOUND."""
    out = np.full((coords.size, 3), np.nan)
    tol = 1.0e-10
    for points, tets, vel in pieces:
        xyz = points[tets]                               # (m, 4, 3)
        lo, hi = xyz.min(axis=1), xyz.max(axis=1)
        keep = np.ones(len(tets), dtype=bool)
        for a, value in fixed.items():
            keep &= (lo[:, a] <= value + tol) & (hi[:, a] >= value - tol)
        if not keep.any():
            continue
        xyz, lo, hi, cand = xyz[keep], lo[keep], hi[keep], tets[keep]
        # BARYCENTRIC MAP: lambda_{1..3} = Tinv (p - v0)
        tmat = np.stack([xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0],
                         xyz[:, 3] - xyz[:, 0]], axis=2)
        tinv = np.linalg.inv(tmat)
        for i, c in enumerate(coords):
            if not np.isnan(out[i, 0]):
                continue
            p = np.empty(3)
            p[axis] = c
            for a, value in fixed.items():
                p[a] = value
            near = np.where((lo[:, axis] <= c + tol) & (hi[:, axis] >= c - tol))[0]
            if near.size == 0:
                continue
            lam = np.einsum("kij,kj->ki", tinv[near], p - xyz[near, 0])
            lam0 = 1.0 - lam.sum(axis=1)
            inside = (lam.min(axis=1) >= -1.0e-9) & (lam0 >= -1.0e-9)
            hit = np.where(inside)[0]
            if hit.size == 0:
                continue
            k = near[hit[0]]
            w = np.concatenate(([lam0[hit[0]]], lam[hit[0]]))
            out[i] = w @ vel[cand[k]]
    return out


def profiles(pieces, z_value, samples):
    allpts = np.vstack([p[0] for p in pieces])
    lo, hi = allpts.min(axis=0), allpts.max(axis=0)
    mid = 0.5 * (lo + hi)
    zc = mid[2] if z_value is None else z_value
    eps = 1.0e-9 * (hi - lo).max()
    ys = np.linspace(lo[1] + eps, hi[1] - eps, samples)
    xs = np.linspace(lo[0] + eps, hi[0] - eps, samples)
    u = sample_line(pieces, {0: mid[0], 2: zc}, 1, ys)[:, 0]
    v = sample_line(pieces, {1: mid[1], 2: zc}, 0, xs)[:, 1]
    print("DOMAIN : x[%g, %g] y[%g, %g] z[%g, %g]; PROFILES AT z = %g"
          % (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], zc))
    for name, arr in (("u(y)", u), ("v(x)", v)):
        missing = np.isnan(arr).sum()
        if missing:
            print("WARNING: %d OF %d %s SAMPLES FOUND NO TETRAHEDRON"
                  % (missing, arr.size, name))
    return ys, u, xs, v, lo, hi


def nice_limits(values, fixed, use_fixed):
    """THE FIXED Ghia Re = 100 RANGE, OR ROUNDED DATA LIMITS FOR OTHER CASES."""
    if use_fixed:
        return fixed
    data = np.concatenate([np.ravel(v) for v in values if v is not None])
    data = data[~np.isnan(data)]
    lo, hi = data.min(), data.max()
    pad = 0.05 * (hi - lo)
    return np.floor((lo - pad) * 20) / 20, np.ceil((hi + pad) * 20) / 20


def load_ref(path):
    data = np.loadtxt(path, delimiter=",", comments="#")
    return data[:, :2]


# ================================
# MAIN
# ================================
def main():
    ap = argparse.ArgumentParser(description="CBS3D lid-driven cavity centreline plots")
    ap.add_argument("result")
    ap.add_argument("--previous", default=None)
    ap.add_argument("--re", type=int, default=100)
    ap.add_argument("--ref-u", default=None)
    ap.add_argument("--ref-v", default=None)
    ap.add_argument("--ref-label", default="Reference")
    ap.add_argument("--label", default="CBS3D")
    ap.add_argument("--z", type=float, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--samples", type=int, default=401)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    family = set_style(plt)

    pieces, used = read_result(args.result)
    ys, u, xs, v, lo, hi = profiles(pieces, args.z, args.samples)
    prefix = args.out or "ldc_re%d" % args.re
    np.savetxt(prefix + "_profiles.csv", np.column_stack([ys, u, xs, v]),
               delimiter=",", header="y,u(x=xc),x,v(y=yc)  from " + used,
               fmt="%.8e")

    # ---- REFERENCE DATA
    if args.ref_u or args.ref_v:
        ref_u = load_ref(args.ref_u) if args.ref_u else None
        ref_v = load_ref(args.ref_v) if args.ref_v else None
        ref_label = args.ref_label
    elif args.re in GHIA:
        ref_u, ref_v = GHIA[args.re]["u"], GHIA[args.re]["v"]
        ref_label = "Ghia et al. (1982), Re = %d" % args.re
    else:
        ref_u = ref_v = None
        ref_label = None
        print("NOTE   : NO BUILT-IN GHIA DATA FOR Re = %d; USE --ref-u/--ref-v" % args.re)

    # ---- AGREEMENT AT THE REFERENCE POINTS
    for name, ref, coord, val in (("u", ref_u, ys, u), ("v", ref_v, xs, v)):
        if ref is None:
            continue
        ok = ~np.isnan(val)
        at = np.interp(ref[:, 0], coord[ok], val[ok])
        print("REF %s  : max |CBS3D - ref| = %.4f at %d points" %
              (name, np.abs(at - ref[:, 1]).max(), len(ref)))

    # ---- CONVERGENCE: CHANGE SINCE AN EARLIER STEP
    if args.previous:
        old, _ = read_result(args.previous)
        _, u0, _, v0, _, _ = profiles(old, args.z, args.samples)
        print("CHANGE : max |du| = %.3e, max |dv| = %.3e since %s"
              % (np.nanmax(np.abs(u - u0)), np.nanmax(np.abs(v - v0)),
                 args.previous))

    # ---- u ALONG THE VERTICAL CENTRELINE
    fig, ax = plt.subplots(figsize=(3.6, 3.3))
    style_axes(ax)
    ax.plot(u, ys, color=LINE_COLOUR, linewidth=1.1, label=args.label, zorder=3)
    if ref_u is not None:
        ax.plot(ref_u[:, 1], ref_u[:, 0], linestyle="none", marker="o",
                markersize=4.2, markerfacecolor="none", markeredgewidth=0.9,
                color=MARK_COLOUR, label=ref_label, zorder=4, clip_on=False)
    ax.set_xlabel(r"$u$")
    ax.set_ylabel(r"$y$")
    ax.set_xlim(*nice_limits([u, None if ref_u is None else ref_u[:, 1]],
                             (-0.3, 1.0), args.re == 100 and not args.ref_u))
    ax.set_ylim(lo[1], hi[1])
    legend_above(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig("%s_u.%s" % (prefix, ext), bbox_inches="tight")
    plt.close(fig)

    # ---- v ALONG THE HORIZONTAL CENTRELINE
    fig, ax = plt.subplots(figsize=(3.6, 3.3))
    style_axes(ax)
    ax.plot(xs, v, color=LINE_COLOUR, linewidth=1.1, label=args.label, zorder=3)
    if ref_v is not None:
        ax.plot(ref_v[:, 0], ref_v[:, 1], linestyle="none", marker="s",
                markersize=4.0, markerfacecolor="none", markeredgewidth=0.9,
                color=MARK_COLOUR, label=ref_label, zorder=4, clip_on=False)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$v$")
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(*nice_limits([v, None if ref_v is None else ref_v[:, 1]],
                             (-0.3, 0.2), args.re == 100 and not args.ref_v))
    legend_above(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig("%s_v.%s" % (prefix, ext), bbox_inches="tight")
    plt.close(fig)

    print("FONT   :", family)
    print("WROTE  : %s_u.png/.pdf, %s_v.png/.pdf, %s_profiles.csv"
          % (prefix, prefix, prefix))


if __name__ == "__main__":
    sys.exit(main())
