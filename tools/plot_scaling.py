#!/usr/bin/env python3
"""
plot_scaling.py - turn scaling_results.csv (from scaling_study.sh) into a
strong-scaling speedup + efficiency figure for CBS3D++_SI on Sunbird.

Usage:
    python plot_scaling.py scaling_results.csv [--out scaling_plot.png]
"""
from __future__ import annotations
import argparse, csv, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?", default="scaling_results.csv")
    ap.add_argument("--out", default="scaling_plot.png")
    args = ap.parse_args()

    threads, secs = [], []
    with open(args.csv) as f:
        for row in csv.DictReader(f):
            threads.append(int(row["threads"]))
            secs.append(float(row["seconds"]))
    if not threads:
        sys.exit("no rows in CSV")

    # sort by thread count, baseline = smallest thread count present
    order = sorted(range(len(threads)), key=lambda i: threads[i])
    threads = [threads[i] for i in order]
    secs = [secs[i] for i in order]
    t1 = secs[0]                      # time at fewest threads (ideally 1)
    base = threads[0]
    speedup = [t1 / s for s in secs]
    ideal = [n / base for n in threads]
    eff = [sp / (n / base) * 100.0 for sp, n in zip(speedup, threads)]

    print(f"{'threads':>8}{'seconds':>12}{'speedup':>10}{'ideal':>8}{'eff %':>8}")
    for n, s, sp, idl, e in zip(threads, secs, speedup, ideal, eff):
        print(f"{n:>8}{s:>12.2f}{sp:>10.2f}{idl:>8.1f}{e:>8.0f}")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))

    ax[0].plot(threads, ideal, "k--", lw=1.3, label="ideal (linear)")
    ax[0].plot(threads, speedup, "o-", color="#1f77b4", lw=2, ms=6, label="CBS3D++_SI")
    ax[0].axvline(20, color="#aaaaaa", ls=":", lw=1)
    ax[0].text(20, ax[0].get_ylim()[1]*0.05, " 1 socket", color="#888", fontsize=8)
    ax[0].set_xlabel("OpenMP threads"); ax[0].set_ylabel("speedup vs %d-thread" % base)
    ax[0].set_title("Strong scaling — speedup"); ax[0].grid(alpha=.4); ax[0].legend()

    ax[1].axhline(100, color="k", ls="--", lw=1.3, label="ideal")
    ax[1].plot(threads, eff, "s-", color="#d62728", lw=2, ms=6, label="parallel efficiency")
    ax[1].axvline(20, color="#aaaaaa", ls=":", lw=1)
    ax[1].set_xlabel("OpenMP threads"); ax[1].set_ylabel("efficiency (%)")
    ax[1].set_title("Parallel efficiency"); ax[1].set_ylim(0, 110)
    ax[1].grid(alpha=.4); ax[1].legend()

    peak = max(speedup)
    fig.suptitle(f"CBS3D++_SI OpenMP scaling on Sunbird — peak speedup {peak:.1f}x "
                 f"at {threads[speedup.index(peak)]} threads", fontweight="bold")
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"\n[ok] wrote {args.out}")


if __name__ == "__main__":
    main()
