#!/usr/bin/env python3
# ================================
# FILE      : plot_tmax_history.py
# PURPOSE   : PUBLICATION FIGURES FROM THE MAXIMUM-TEMPERATURE HISTORY OF A
#             PULSED-HEATING RUN (CSV: time_s,iteration,dt_s,Tmax_K).
#
#             FIGURE 1  HEATING FACTOR f(t) AND T_max(t), WITH THE PULSE
#                       PHASES, THE STARTING STEADY STATE, THE FINAL PLATEAU
#                       AND THE COLDEST POINT OF THE DWELL MARKED.
#             FIGURE 2  NORMALISED COOLING AND REHEATING RESPONSES WITH
#                       SINGLE-EXPONENTIAL FITS (TIME CONSTANTS).
#
#             THE FITS T(t) = T_inf + A exp(-(t - t_s)/tau) USE THE DATA FROM
#             10 s AFTER THE POWER REACHES ITS NEW LEVEL (t_s) TO THE END OF
#             THE PHASE, SO THE SHORT LAG RIGHT AFTER THE SWITCH IS EXCLUDED.
#
# USAGE     : python3 tools/plot_tmax_history.py Tmax_history.csv
#                 [--pulse cases/blanket_slice/pulse_600s.txt]
#                 [--t-steady 437.06] [--t-inlet 300] [--out figures/tmax]
#
# OUTPUT    : <out>_history.pdf/.png, <out>_time_constants.pdf/.png,
#             <out>_summary.txt (NUMBERS QUOTED IN THE FIGURES)
# ================================

from __future__ import print_function

import argparse
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

# ================================
# STYLE: WHITE BACKGROUND, 9 pt TEXT, THIN RECESSIVE AXES
# ================================
INK = "#1a1a1a"
INK2 = "#55554f"
GRID = "#e4e4e0"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
DWELL = "#ececea"
RAMP = "#f6f6f4"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.edgecolor": INK2,
    "axes.labelcolor": INK,
    "axes.linewidth": 0.8,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelcolor": INK,
    "ytick.labelcolor": INK,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "legend.fontsize": 8.5,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "pdf.fonttype": 42,
})


def read_history(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    t = np.asarray(d["time_s"], float)
    T = np.asarray(d["Tmax_K"], float)
    o = np.argsort(t, kind="stable")
    t, T = t[o], T[o]
    keep = np.r_[True, np.diff(t) > 0.0]
    return t[keep], T[keep]


def read_pulse(path):
    rows = []
    for line in open(path):
        s = line.split("#")[0].strip()
        if s:
            rows.append([float(v) for v in s.split()[:2]])
    return np.array(rows)


def fit_phase(t, T, t_switch, t_end, lag=10.0, p0=None):
    m = (t >= t_switch + lag) & (t <= t_end)

    def model(tt, Tinf, A, tau):
        return Tinf + A * np.exp(-(tt - t_switch) / tau)

    p, _ = curve_fit(model, t[m], T[m], p0=p0, maxfev=20000)
    r = T[m] - model(t[m], *p)
    return dict(Tinf=p[0], A=p[1], tau=p[2], rms=float(np.sqrt(np.mean(r ** 2))),
                maxres=float(np.abs(r).max()), model=model, p=p,
                t_switch=t_switch, t_end=t_end, lag=lag)


def shade_phases(ax, pulse):
    # RAMP DOWN, DWELL, RAMP UP FROM THE PULSE TABLE (f = 1 -> 0 -> 1)
    tt, ff = pulse[:, 0], pulse[:, 1]
    for i in range(len(tt) - 1):
        a, b = tt[i], tt[i + 1]
        if ff[i] == 0.0 and ff[i + 1] == 0.0:
            ax.axvspan(a, b, color=DWELL, lw=0, zorder=0)
        elif ff[i] != ff[i + 1]:
            ax.axvspan(a, b, color=RAMP, lw=0, zorder=0)


def style_axes(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis="y", color=GRID, lw=0.6, zorder=0.5)
    ax.set_axisbelow(True)


def figure_history(t, T, pulse, cool, heat, T_steady, T_inlet, out):
    t_end = pulse[-1, 0]
    i_min = int(np.argmin(T))
    plateau = heat["Tinf"]

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(7.0, 4.9), sharex=True,
        gridspec_kw=dict(height_ratios=[1, 3.2], hspace=0.08))

    # ---- (a) HEATING FACTOR
    shade_phases(ax0, pulse)
    ax0.plot(pulse[:, 0], pulse[:, 1], color=INK, lw=1.5, solid_joinstyle="miter")
    ax0.set_ylim(-0.12, 1.25)
    ax0.set_yticks([0, 1])
    ax0.set_yticklabels(["0", "1"])
    ax0.set_ylabel("heating\nfactor $f$")
    style_axes(ax0)
    ax0.grid(False)
    ax0.text(5, 1.08, "full power", ha="left", va="bottom", fontsize=8, color=INK2)
    ax0.text(150, 0.12, "dwell: no heating", ha="center", va="bottom", fontsize=8, color=INK2)
    ax0.text(435, 1.08, "full power", ha="center", va="bottom", fontsize=8, color=INK2)
    ax0.text(-0.115, 1.0, "(a)", transform=ax0.transAxes, fontsize=10, fontweight="bold", va="top")

    # ---- (b) MAXIMUM TEMPERATURE
    shade_phases(ax1, pulse)
    style_axes(ax1)
    ax1.axhline(T_steady, color=INK2, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax1.axhline(plateau, color=INK2, lw=0.9, ls=(0, (1, 2)), zorder=1)
    ax1.axhline(T_inlet, color=INK2, lw=0.9, ls=(0, (6, 2, 1, 2)), zorder=1)
    ax1.plot(t, T, color=BLUE, lw=1.6, zorder=3, solid_joinstyle="round")

    ax1.text(44, T_steady + 1.5, "start: steady state of run 30102, %.0f K" % T_steady,
             ha="left", va="bottom", fontsize=8, color=INK2)
    ax1.text(t_end - 4, plateau + 2.0, "pulsed-run plateau, %.0f K" % plateau,
             ha="right", va="bottom", fontsize=8, color=INK2)
    ax1.text(t_end - 4, T_inlet + 1.5, "coolant inlet, %.0f K" % T_inlet,
             ha="right", va="bottom", fontsize=8, color=INK2)

    # COLDEST POINT OF THE DWELL
    ax1.plot(t[i_min], T[i_min], "o", ms=5.5, mfc="white", mec=BLUE, mew=1.4, zorder=4)
    ax1.annotate("minimum %.0f K at t = %.0f s" % (T[i_min], t[i_min]),
                 xy=(t[i_min], T[i_min]), xytext=(t[i_min] + 52, T[i_min] + 10),
                 ha="left", va="center", fontsize=8, color=INK,
                 arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7,
                                 shrinkA=2, shrinkB=4))

    # TIME CONSTANTS ON THE CURVE
    ax1.text(102, 372, r"cooling, $\tau \approx %.0f$ s" % cool["tau"],
             ha="left", va="center", fontsize=8, color=INK)
    ax1.text(318, 372, r"reheating, $\tau \approx %.0f$ s" % heat["tau"],
             ha="left", va="center", fontsize=8, color=INK)

    # OFFSET BETWEEN THE START AND THE PLATEAU
    x_arrow = 588.0
    ax1.annotate("", xy=(x_arrow, plateau), xytext=(x_arrow, T_steady),
                 arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9,
                                 shrinkA=0, shrinkB=0))
    ax1.text(x_arrow - 5, 0.5 * (plateau + T_steady), "+%.0f K" % (plateau - T_steady),
             ha="right", va="center", fontsize=8.5, color=INK, fontweight="bold")

    ax1.set_xlim(0, t_end)
    ax1.set_ylim(290, 462)
    ax1.set_xticks(np.arange(0, t_end + 1, 60))
    ax1.set_xlabel("time  $t$  [s]")
    ax1.set_ylabel(r"maximum temperature,  $T_{\max}$  [K]")
    ax1.text(-0.115, 1.0, "(b)", transform=ax1.transAxes, fontsize=10, fontweight="bold", va="top")

    for ext in ("pdf", "png"):
        fig.savefig(out + "_history." + ext)
    plt.close(fig)


def figure_time_constants(t, T, cool, heat, out):
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    style_axes(ax)
    ax.grid(True, axis="both", color=GRID, lw=0.6)

    e = np.exp(-1.0)
    rows = [(cool, BLUE, "cooling: heating off at t = %.0f s" % cool["t_switch"]),
            (heat, ORANGE, "reheating: full power from t = %.0f s" % heat["t_switch"])]
    for ph, col, lab in rows:
        m = (t >= ph["t_switch"]) & (t <= ph["t_end"])
        s = t[m] - ph["t_switch"]
        theta = (T[m] - ph["Tinf"]) / ph["A"]
        ax.plot(s, theta, color=col, lw=1.6, zorder=3, label=lab)
        sf = np.linspace(ph["lag"], ph["t_end"] - ph["t_switch"], 200)
        thf = np.exp(-sf / ph["tau"])
        ax.plot(sf, thf, color=INK, lw=0.9, ls=(0, (4, 2.5)), zorder=4)
        # THE TIME CONSTANT: WHERE THE FIT HAS FALLEN BY A FACTOR e
        th_tau = e
        ax.plot([ph["tau"], ph["tau"]], [-0.04, th_tau], color=col, lw=0.9, ls=(0, (1, 1.5)), zorder=2)
        ax.plot(ph["tau"], th_tau, "o", ms=5, mfc="white", mec=col, mew=1.3, zorder=5)

    ax.plot([], [], color=INK, lw=0.9, ls=(0, (4, 2.5)), label=r"fit, $\exp(-s/\tau)$")
    ax.text(cool["tau"] - 3, 0.03, r"$\tau$ = %.0f s" % cool["tau"], ha="right", va="bottom",
            fontsize=8, color=INK)
    ax.text(heat["tau"] + 3, 0.03, r"$\tau$ = %.0f s" % heat["tau"], ha="left", va="bottom",
            fontsize=8, color=INK)
    ax.axhline(e, color=INK2, lw=0.6, ls=(0, (1, 2)), zorder=1)
    ax.text(318, e + 0.015, "1/e", ha="right", va="bottom", fontsize=8, color=INK2)

    ax.set_xlim(0, 330)
    ax.set_ylim(-0.04, 1.12)
    ax.set_xticks(np.arange(0, 331, 30))
    ax.set_xlabel("time since the heating reached its new level,  $s$  [s]")
    ax.set_ylabel("change still to come,  "
                  r"$(T_{\max}-T_\infty)\,/\,\Delta T$")
    ax.legend(loc="upper right", handlelength=2.2)

    for ext in ("pdf", "png"):
        fig.savefig(out + "_time_constants." + ext)
    plt.close(fig)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--pulse", default=os.path.join(here, "..", "cases", "blanket_slice", "pulse_600s.txt"))
    ap.add_argument("--t-steady", type=float, default=None,
                    help="T_max OF THE STARTING STEADY STATE [K] (DEFAULT: FIRST SAMPLE)")
    ap.add_argument("--t-inlet", type=float, default=300.0)
    ap.add_argument("--out", default="tmax")
    args = ap.parse_args()

    t, T = read_history(args.csv)
    pulse = read_pulse(args.pulse)
    T_steady = args.t_steady if args.t_steady is not None else T[0]

    # PHASE TIMES FROM THE PULSE TABLE: HEATING REACHES 0, THEN REACHES 1 AGAIN
    tp, fp = pulse[:, 0], pulse[:, 1]
    i0 = int(np.where(fp == 0.0)[0][0])
    i1 = int(np.where((fp == 1.0) & (np.arange(len(fp)) > i0))[0][0])
    t_off, t_dwell_end, t_on = tp[i0], tp[i1 - 1], tp[i1]

    cool = fit_phase(t, T, t_off, t_dwell_end, p0=(args.t_inlet, 90.0, 50.0))
    heat = fit_phase(t, T, t_on, t[-1], p0=(T.max(), -100.0, 60.0))

    d = os.path.dirname(args.out)
    if d:
        os.makedirs(d, exist_ok=True)
    figure_history(t, T, pulse, cool, heat, T_steady, args.t_inlet, args.out)
    figure_time_constants(t, T, cool, heat, args.out)

    i_min = int(np.argmin(T))
    with open(args.out + "_summary.txt", "w") as f:
        f.write("samples                 %d  (t = %.3f to %.3f s)\n" % (len(t), t[0], t[-1]))
        f.write("T_max at t = 0          %.1f K\n" % T[0])
        f.write("T_max, starting steady  %.2f K\n" % T_steady)
        f.write("T_max minimum           %.1f K at t = %.1f s\n" % (T[i_min], t[i_min]))
        f.write("T_max, last 60 s        mean %.2f K, range %.0f-%.0f K\n" % (
            T[t >= t[-1] - 60].mean(), T[t >= t[-1] - 60].min(), T[t >= t[-1] - 60].max()))
        for name, ph in (("cooling", cool), ("reheating", heat)):
            f.write("%-10s fit window %.0f-%.0f s: tau = %.1f s, T_inf = %.2f K, "
                    "rms residual %.2f K, max %.2f K\n" % (
                        name, ph["t_switch"] + ph["lag"], ph["t_end"], ph["tau"],
                        ph["Tinf"], ph["rms"], ph["maxres"]))
    print(open(args.out + "_summary.txt").read())


if __name__ == "__main__":
    main()
