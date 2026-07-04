#!/usr/bin/env python3
"""
CBS3D++_SI live residual plotter.

Usage:
    python tools/live_residual_plot.py output/<case>/<case>_residuals.csv

The solver writes a CSV file.  This script watches it and updates one semilog
figure containing all important residual histories.
"""

from __future__ import annotations

import argparse
import csv
import math
import time
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt


SERIES = [
    ("u_rel", "u"),
    ("v_rel", "v"),
    ("w_rel", "w"),
    ("p_rel", "p"),
    ("T_rel", "T"),
    ("velocity_rel_max", "max velocity"),
    ("cg_relative_l2", "CG relative L2"),
]


def read_csv(path: Path) -> Dict[str, List[float]]:
    if not path.exists():
        return {}

    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        return {}

    data: Dict[str, List[float]] = {key: [] for key in rows[0].keys()}

    for row in rows:
        for key, value in row.items():
            try:
                data[key].append(float(value))
            except (TypeError, ValueError):
                data[key].append(float("nan"))

    return data


def positive(values: List[float]) -> List[float]:
    cleaned: List[float] = []
    tiny = 1.0e-300

    for value in values:
        if math.isfinite(value) and value > 0.0:
            cleaned.append(value)
        else:
            cleaned.append(tiny)

    return cleaned


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--interval", type=float, default=1.0)
    args = parser.parse_args()

    plt.ion()
    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    fig.canvas.manager.set_window_title("CBS3D++_SI residual monitor")

    while True:
        data = read_csv(args.csv_file)

        ax.clear()
        ax.set_title("CBS3D++_SI residual history")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Residual")
        ax.set_yscale("log")
        ax.grid(True, which="major", linestyle="-", linewidth=0.6, alpha=0.65)
        ax.grid(True, which="minor", linestyle=":", linewidth=0.45, alpha=0.55)
        ax.minorticks_on()

        if data and "iteration" in data:
            iterations = data["iteration"]

            for key, label in SERIES:
                if key in data:
                    ax.plot(iterations, positive(data[key]), label=label, linewidth=1.8)

            ax.legend(loc="best", frameon=True)
            ax.set_xlim(left=0)
        else:
            ax.text(
                0.5,
                0.5,
                f"Waiting for residual CSV:\n{args.csv_file}",
                transform=ax.transAxes,
                ha="center",
                va="center",
            )

        fig.tight_layout()
        plt.pause(args.interval)


if __name__ == "__main__":
    main()
