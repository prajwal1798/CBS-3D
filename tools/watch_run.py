#!/usr/bin/env python3
# ================================
# FILE      : watch_run.py
# PURPOSE   : LIVE TERMINAL MONITOR OF A CBS3D MPI RUN (STEADY OR TRANSIENT).
#             REFRESHES IN PLACE EVERY --every SECONDS FROM
#               output/<case>/<case>_distributed_residuals.csv  (EVERY 10 ITS)
#               solver.out                                      (STOP REASON, ERRORS)
#               squeue                                          (JOB STATE)
#
#             SHOWS: ITERATION, SPEED, ELAPSED; EACH RELATIVE RESIDUAL WITH
#             ITS DECADES PER 1000 ITERATIONS AND STATUS AGAINST THE STOP
#             TARGET; A HISTORY OF log10(max velocity residual); AN ESTIMATE
#             OF THE ITERATIONS AND TIME TO THE TARGET (LOG-LINEAR FIT OVER
#             THE LAST QUARTER OF THE HISTORY); CONTINUITY, CG, max|u|, SA.
#
# USAGE     : python3 tools/watch_run.py RUN_DIR [--every 5] [--u-target 1e-7]
#                 [--sa-target 1e-6] [--unicode] [--no-color] [--once]
#             (--unicode DRAWS THE HISTORY WITH BLOCK CHARACTERS; PuTTY'S
#             DEFAULT FONT SHOWS THEM AS BOXES, SO ASCII IS THE DEFAULT)
#             RUN_DIR = THE run_<case>_<job> FOLDER. Ctrl+C STOPS WATCHING
#             (THE JOB CONTINUES).
# ================================

from __future__ import print_function

import argparse
import glob
import math
import os
import re
import subprocess
import sys
import time

BLOCKS_UNI = "▁▂▃▄▅▆▇█"
BLOCKS_ASC = "_.,-~=+*#"


class Csv(object):
    # INCREMENTAL READER: ONLY NEW COMPLETE LINES ARE PARSED ON EACH CALL
    def __init__(self, path):
        self.path = path
        self.pos = 0
        self.head = None
        self.rows = []

    def update(self):
        if not os.path.exists(self.path):
            return
        with open(self.path, "r") as f:
            f.seek(self.pos)
            chunk = f.read()
        if not chunk:
            return
        last = chunk.rfind("\n")
        if last < 0:
            return
        self.pos += len(chunk[:last + 1].encode())
        for line in chunk[:last].split("\n"):
            if not line.strip():
                continue
            if self.head is None:
                self.head = line.strip().split(",")
                continue
            try:
                self.rows.append([float(v) for v in line.split(",")])
            except ValueError:
                pass

    def col(self, name):
        if self.head is None or name not in self.head:
            return None
        i = self.head.index(name)
        return [r[i] for r in self.rows if len(r) > i]


def fmt_time(s):
    if s is None or not math.isfinite(s) or s < 0:
        return "--:--:--"
    s = int(round(s))
    return "%d:%02d:%02d" % (s // 3600, (s % 3600) // 60, s % 60)


def fit_slope(its, vals, frac=0.25, nmin=20):
    # SLOPE OF log10(val) AGAINST ITERATION OVER THE LAST frac OF THE HISTORY
    pts = [(i, math.log10(v)) for i, v in zip(its, vals) if v > 0 and math.isfinite(v)]
    if len(pts) < nmin:
        return None
    pts = pts[-max(nmin, int(len(pts) * frac)):]
    n = float(len(pts))
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    if sxx <= 0:
        return None
    return sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx


def sparkline(vals, width, blocks):
    vals = [math.log10(v) for v in vals if v > 0 and math.isfinite(v)]
    if not vals:
        return ""
    step = max(1, len(vals) // width)
    bins = [min(vals[i:i + step]) for i in range(0, len(vals), step)][-width:]
    lo, hi = min(bins), max(bins)
    span = (hi - lo) or 1.0
    line = "".join(blocks[int((b - lo) / span * (len(blocks) - 1))] for b in bins)
    return line, lo, hi


def job_state(run_dir):
    m = re.search(r"_(\d+)$", os.path.basename(os.path.normpath(run_dir)))
    if not m:
        return None, "?"
    jid = m.group(1)
    try:
        out = subprocess.run(["squeue", "-h", "-j", jid, "-o", "%T %M"], capture_output=True,
                             text=True, timeout=5).stdout.strip()
        return jid, out if out else "NOT IN QUEUE (FINISHED)"
    except Exception:
        return jid, "squeue unavailable"


def log_tail(run_dir):
    path = os.path.join(run_dir, "solver.out")
    if not os.path.exists(path):
        return None, []
    with open(path, "rb") as f:
        f.seek(0, 2)
        f.seek(max(0, f.tell() - 200000))
        text = f.read().decode("utf-8", "replace")
    text = re.sub(r"\x1b\[[0-9;]*m", "", text).replace("\r", "\n")
    stop = None
    m = re.search(r"stop reason\s*:\s*(.+)", text)
    if m:
        stop = m.group(1).strip()
    errs = [ln.strip() for ln in text.split("\n") if re.search(r"ERROR|FATAL|NaN|nan\b", ln)][-3:]
    return stop, errs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--every", type=float, default=5.0)
    ap.add_argument("--u-target", type=float, default=1.0e-7)
    ap.add_argument("--sa-target", type=float, default=1.0e-6)
    ap.add_argument("--unicode", action="store_true")
    ap.add_argument("--no-color", action="store_true")
    ap.add_argument("--once", action="store_true")
    args = ap.parse_args()

    run = os.path.abspath(args.run_dir)
    found = sorted(glob.glob(os.path.join(run, "output", "*", "*_residuals.csv")))
    while not found and not args.once:
        print("WAITING FOR %s/output/*/*_residuals.csv ..." % run)
        time.sleep(args.every)
        found = sorted(glob.glob(os.path.join(run, "output", "*", "*_residuals.csv")))
    if not found:
        sys.exit("NO RESIDUAL FILE UNDER %s/output" % run)

    csv = Csv(found[0])
    blocks = BLOCKS_UNI if args.unicode else BLOCKS_ASC
    color = sys.stdout.isatty() and not args.no_color
    G, Y, R, B, D, N = (("\033[32m", "\033[33m", "\033[31m", "\033[1m", "\033[2m", "\033[0m")
                        if color else ("",) * 6)
    t_start = time.time()

    while True:
        csv.update()
        its = csv.col("iteration") or []
        lines = []
        bar = "=" * 74
        jid, state = job_state(run)
        lines.append(B + "CBS3D RUN MONITOR" + N + "   %s   job %s: %s   %s" % (
            os.path.basename(run), jid, state, time.strftime("%H:%M:%S")))
        lines.append(bar)

        if its:
            wall = csv.col("iteration_wall_seconds") or []
            recent = wall[-50:]
            s_per_it = sum(recent) / len(recent) if recent else float("nan")
            stride = its[1] - its[0] if len(its) > 1 else 1.0     # LOGGED EVERY stride ITERATIONS
            elapsed = sum(wall) * stride
            lines.append("iteration  %s%10d%s     %6.1f it/s     elapsed %s     dt %s" % (
                B, int(its[-1]), N, 1.0 / s_per_it if s_per_it > 0 else 0.0, fmt_time(elapsed),
                "%.3e" % csv.col("dt")[-1] if csv.col("dt") else "-"))
            lines.append("")
            lines.append("%-12s %12s %16s %12s   %s" % ("residual", "current", "dec / 1000 it", "target", "status"))
            lines.append("-" * 74)

            def row(label, name, target):
                v = csv.col(name)
                if not v:
                    return
                s = fit_slope(its, v)
                rate = "%+.2f" % (1000.0 * s) if s is not None else "-"
                if target is None:
                    st, c = "", D
                elif v[-1] <= target * 1.005:          # PRINT ROUNDING AT THE STOP
                    st, c = "MET", G
                else:
                    st, c = "%.1f dec to go" % math.log10(v[-1] / target), Y
                lines.append("%-12s %s%12.3e%s %16s %12s   %s%s%s" % (
                    label, c, v[-1], N, rate, "%.0e" % target if target else "", c, st, N))

            row("velocity", "velocity_rel_max", args.u_target)
            row("  u", "u_rel", None)
            row("  v", "v_rel", None)
            row("  w", "w_rel", None)
            row("pressure", "p_rel", None)
            if csv.col("sa_rel"):
                row("SA nu~", "sa_rel", args.sa_target)
            if (csv.col("T_rel") or [0])[-1] > 0:
                row("temperature", "T_rel", None)

            vel = csv.col("velocity_rel_max")
            lines.append("")
            sp = sparkline(vel, 60, blocks)
            if sp:
                lines.append("history log10(velocity residual)  [%+.1f .. %+.1f]" % (sp[2], sp[1]))
                lines.append("  " + sp[0])
                lines.append("  it %-10d%s%10d" % (int(its[0]), " " * 38, int(its[-1])))

            s = fit_slope(its, vel)
            if vel and vel[-1] > args.u_target and s is not None and s < 0:
                need = math.log10(args.u_target / vel[-1]) / s
                lines.append("estimate   %d more iterations, about %s (log-linear fit, last quarter)" % (
                    int(need), fmt_time(need * s_per_it)))
            elif vel and vel[-1] > args.u_target:
                lines.append("estimate   " + Y + "not converging at present (residual flat or rising)" + N)
            lines.append("")

            def last(name, f="%.3e"):
                v = csv.col(name)
                return f % v[-1] if v else "-"
            lines.append("continuity rms %s   max %s        max|u| %s   max du %s" % (
                last("continuity_rms"), last("continuity_max"), last("maximum_velocity", "%.4f"),
                last("maximum_velocity_correction")))
            lines.append("pressure   CG %s its, rel %s" % (last("cg_iterations", "%d"), last("cg_relative_l2")))
            if csv.col("sa_rel"):
                lines.append("SA         nu~ [%s, %s]   max nu_t %s   max chi %s" % (
                    last("sa_nu_tilde_min"), last("sa_nu_tilde_max"), last("sa_nu_t_max"),
                    last("sa_chi_max", "%.1f")))
        else:
            lines.append("NO ITERATIONS YET (%s)" % csv.path)

        stop, errs = log_tail(run)
        lines.append(bar)
        if stop:
            lines.append(G + B + "FINISHED: " + stop + N)
        for e in errs:
            lines.append(R + e[:120] + N)
        lines.append(D + "refresh %gs   Ctrl+C stops watching (the job continues)" % args.every + N)

        if args.once:
            print("\n".join(lines))
            return
        sys.stdout.write("\033[H\033[2J" + "\n".join(lines) + "\n")
        sys.stdout.flush()
        if stop:
            return
        try:
            time.sleep(args.every)
        except KeyboardInterrupt:
            print()
            return


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print()
