#!/usr/bin/env python3
# ================================
# FILE      : plot_residuals_pub.py
# PURPOSE   : CBS3D CONVERGENCE HISTORY (u, v, w, p AND, IF SOLVED, T RELATIVE
#             CHANGES PER ITERATION) IN THE SAME STYLE AS
#             plot_ldc_centerlines.py, PLUS A PLAIN-WORDS SUMMARY:
#             PSEUDO-TIME REACHED, TIME STEP, AND HOW FAST THE RESIDUALS
#             ARE STILL FALLING.
#
# USAGE     : python plot_residuals_pub.py RESIDUALS [--out PREFIX]
#
#   RESIDUALS <case>_residuals.csv, <case>_distributed_residuals.csv, OR THE
#             output/<case> FOLDER THAT HOLDS ONE OF THEM.
#
# OUTPUT    : PREFIX.png (600 dpi) AND PREFIX.pdf   (DEFAULT PREFIX: residuals)
#
# NOTE      : u_rel IS THE CHANGE OF u IN ONE ITERATION. WITH A SMALL TIME STEP
#             IT IS SMALL EVEN WHILE THE FLOW IS STILL DEVELOPING; u_rel / dt
#             (A RATE) AND THE PSEUDO-TIME REACHED ARE THE BETTER GUIDES.
# ================================

from __future__ import print_function

import argparse
import csv
import os
import sys

import numpy as np

BACKGROUND = "#F3F0FA"
FONTS = ["Arial", "Liberation Sans", "Helvetica", "Arimo", "DejaVu Sans"]
SERIES = [("u_rel", r"$u$", "#1F4E9C"), ("v_rel", r"$v$", "#E0661B"),
          ("w_rel", r"$w$", "#2A9D5B"), ("p_rel", r"$p$", "#8E3BA8"),
          ("T_rel", r"$T$", "#C0392B")]


def set_style(plt):
    from matplotlib import font_manager
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
        "font.family": "sans-serif", "font.sans-serif": [family],
        "mathtext.fontset": "custom", "mathtext.rm": family,
        "mathtext.it": family + ":italic", "mathtext.bf": family + ":bold",
        "font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8.5,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.linewidth": 0.8, "axes.edgecolor": "black",
        "xtick.color": "black", "ytick.color": "black",
        "xtick.major.size": 3.5, "ytick.major.size": 3.5,
        "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
        "savefig.dpi": 600, "pdf.fonttype": 42,
    })
    return family


def find_csv(path):
    if os.path.isdir(path):
        files = sorted(f for f in os.listdir(path) if f.endswith("_residuals.csv"))
        if not files:
            raise ValueError("NO *_residuals.csv IN %s" % path)
        path = os.path.join(path, files[-1])
    return path


def read_csv(path):
    with open(path) as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("%s HAS NO DATA ROWS" % path)
    cols = {}
    for key in rows[0]:
        try:
            cols[key] = np.array([float(r[key]) for r in rows])
        except (TypeError, ValueError):
            pass
    return cols


def main():
    ap = argparse.ArgumentParser(description="CBS3D residual history")
    ap.add_argument("residuals")
    ap.add_argument("--out", default="residuals")
    args = ap.parse_args()

    path = find_csv(args.residuals)
    c = read_csv(path)
    it, dt = c["iteration"], c["dt"]

    # PSEUDO-TIME: EACH LOGGED ROW STANDS FOR THE ITERATIONS SINCE THE LAST ONE.
    steps = np.diff(np.concatenate(([0.0], it)))
    ptime = np.cumsum(steps * dt)

    print("READ   :", path)
    print("RUN    : %d iterations, dt now %.3e (min %.3e, max %.3e)"
          % (it[-1], dt[-1], dt.min(), dt.max()))
    print("TIME   : pseudo-time reached = %.3f (lid-velocity / cavity-size units)"
          % ptime[-1])

    # HOW FAST ARE THE RESIDUALS STILL FALLING? COMPARE THE LAST 10 % OF THE RUN
    # WITH THE 10 % BEFORE IT.
    n = len(it)
    k = max(1, n // 10)
    if n >= 3 * k:
        for key, name, _ in SERIES[:4]:
            if key in c:
                late = np.median(c[key][-k:])
                early = np.median(c[key][-2 * k:-k])
                span = it[-1] - it[-2 * k]
                print("TREND  : %s_rel %.2e -> %.2e over the last ~%d iterations (x%.2f)"
                      % (key[0], early, late, span, late / early if early else np.nan))
        print("RATE   : u_rel/dt = %.2e   (the flow is still changing at this rate per unit time)"
              % (c["u_rel"][-1] / dt[-1]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter
    family = set_style(plt)

    fig, ax = plt.subplots(figsize=(4.4, 3.3))
    ax.set_facecolor(BACKGROUND)
    ax.set_axisbelow(True)
    for key, name, colour in SERIES:
        if key in c and np.any(c[key] > 0):
            ax.semilogy(it, c[key], color=colour, linewidth=1.0, label=name)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Relative change per iteration")
    ax.set_xlim(0, it[-1])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=15))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=15))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(which="major", color="white", linewidth=0.9)
    ax.grid(which="minor", color="white", linewidth=0.45)
    ax.tick_params(which="both", top=True, right=True)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=5,
              frameon=False, handlelength=2.0, columnspacing=1.4)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig("%s.%s" % (args.out, ext), bbox_inches="tight")
    print("FONT   :", family)
    print("WROTE  : %s.png, %s.pdf" % (args.out, args.out))


if __name__ == "__main__":
    sys.exit(main())
