#!/usr/bin/env python3
"""Analyze CBS3D MPI strong-scaling runs and produce publication-ready metrics.

The smallest available rank count is used as the measured baseline. For the
initial P40/P80 campaign this means that speedup and efficiency are relative to
P40, not to one core. The same script extends directly to P160--P1280.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path
from typing import Dict, List


def read_key_values(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if "=" in raw:
            key, value = raw.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def latest_complete_job(level_dir: Path) -> Path:
    candidates = []
    for job_dir in level_dir.glob("job_*"):
        if (job_dir / "timings.csv").is_file() and (job_dir / "statistics.txt").is_file():
            try:
                job_id = int(job_dir.name.split("_", 1)[1])
            except (IndexError, ValueError):
                continue
            candidates.append((job_id, job_dir))
    if not candidates:
        raise FileNotFoundError(f"no completed benchmark job under {level_dir}")
    return max(candidates, key=lambda item: item[0])[1]


def load_level(level_dir: Path) -> Dict[str, float]:
    job_dir = latest_complete_job(level_dir)
    metadata = read_key_values(job_dir / "metadata.txt")

    with (job_dir / "timings.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    measured = [row for row in rows if row["classification"] == "measured"]
    if not measured:
        raise ValueError(f"no measured rows in {job_dir / 'timings.csv'}")

    loop_times = [float(row["loop_time_s"]) for row in measured]
    setup_times = [float(row["setup_time_s"]) for row in measured]
    elapsed_times = [float(row["elapsed_time_s"]) for row in measured]

    ranks = int(metadata["mpi_ranks"])
    nodes = int(metadata["nodes"])
    iterations = int(metadata["iterations"])

    return {
        "ranks": ranks,
        "nodes": nodes,
        "iterations": iterations,
        "measured_runs": len(measured),
        "median_loop_time_s": statistics.median(loop_times),
        "minimum_loop_time_s": min(loop_times),
        "maximum_loop_time_s": max(loop_times),
        "mean_loop_time_s": statistics.mean(loop_times),
        "stddev_loop_time_s": statistics.pstdev(loop_times),
        "median_setup_time_s": statistics.median(setup_times),
        "median_elapsed_time_s": statistics.median(elapsed_times),
        "job_id": int(job_dir.name.split("_", 1)[1]),
        "job_dir": str(job_dir),
    }


def write_plots(rows: List[Dict[str, float]], plot_dir: Path) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    plot_dir.mkdir(parents=True, exist_ok=True)
    ranks = [int(row["ranks"]) for row in rows]

    def save_current(name: str) -> None:
        plt.tight_layout()
        plt.savefig(plot_dir / f"{name}.png", dpi=300, bbox_inches="tight")
        plt.savefig(plot_dir / f"{name}.pdf", bbox_inches="tight")
        plt.close()

    medians = [row["median_loop_time_s"] for row in rows]
    lower = [row["median_loop_time_s"] - row["minimum_loop_time_s"] for row in rows]
    upper = [row["maximum_loop_time_s"] - row["median_loop_time_s"] for row in rows]
    plt.figure(figsize=(6.4, 4.4))
    plt.errorbar(ranks, medians, yerr=[lower, upper], marker="o", capsize=4)
    plt.xlabel("MPI ranks")
    plt.ylabel("Median solver-loop time (s)")
    plt.title("CBS3D strong scaling: runtime")
    plt.grid(True, alpha=0.3)
    save_current("runtime_vs_ranks")

    speedups = [row["relative_speedup"] for row in rows]
    ideals = [row["ideal_relative_speedup"] for row in rows]
    plt.figure(figsize=(6.4, 4.4))
    plt.plot(ranks, speedups, marker="o", label="Measured")
    plt.plot(ranks, ideals, linestyle="--", label="Ideal")
    plt.xlabel("MPI ranks")
    plt.ylabel("Speedup relative to baseline")
    plt.title("CBS3D strong scaling: speedup")
    plt.grid(True, alpha=0.3)
    plt.legend()
    save_current("speedup_vs_ranks")

    efficiencies = [row["relative_parallel_efficiency_pct"] for row in rows]
    plt.figure(figsize=(6.4, 4.4))
    plt.plot(ranks, efficiencies, marker="o")
    plt.axhline(100.0, linestyle="--")
    plt.xlabel("MPI ranks")
    plt.ylabel("Parallel efficiency relative to baseline (%)")
    plt.title("CBS3D strong scaling: parallel efficiency")
    plt.grid(True, alpha=0.3)
    save_current("parallel_efficiency_vs_ranks")

    throughput = [row["iterations_per_second"] for row in rows]
    plt.figure(figsize=(6.4, 4.4))
    plt.plot(ranks, throughput, marker="o")
    plt.xlabel("MPI ranks")
    plt.ylabel("Iterations per second")
    plt.title("CBS3D strong scaling: throughput")
    plt.grid(True, alpha=0.3)
    save_current("throughput_vs_ranks")

    core_hours = [row["core_hours_per_solve"] for row in rows]
    plt.figure(figsize=(6.4, 4.4))
    plt.plot(ranks, core_hours, marker="o")
    plt.xlabel("MPI ranks")
    plt.ylabel("Core-hours per 10,000-iteration solve")
    plt.title("CBS3D strong scaling: computational cost")
    plt.grid(True, alpha=0.3)
    save_current("core_hours_vs_ranks")

    setup_fraction = [row["setup_fraction_pct"] for row in rows]
    plt.figure(figsize=(6.4, 4.4))
    plt.plot(ranks, setup_fraction, marker="o")
    plt.xlabel("MPI ranks")
    plt.ylabel("Setup fraction of setup + loop time (%)")
    plt.title("CBS3D strong scaling: setup overhead")
    plt.grid(True, alpha=0.3)
    save_current("setup_fraction_vs_ranks")

    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-root", type=Path, required=True)
    args = parser.parse_args()

    level_dirs = sorted(path for path in args.campaign_root.glob("P[0-9][0-9][0-9][0-9]") if path.is_dir())
    if not level_dirs:
        raise SystemExit(f"no partition-level results under {args.campaign_root}")

    rows = [load_level(level_dir) for level_dir in level_dirs]
    rows.sort(key=lambda row: int(row["ranks"]))

    baseline_ranks = int(rows[0]["ranks"])
    baseline_time = float(rows[0]["median_loop_time_s"])

    for row in rows:
        ranks = int(row["ranks"])
        median_time = float(row["median_loop_time_s"])
        iterations = int(row["iterations"])
        ideal_speedup = ranks / baseline_ranks
        speedup = baseline_time / median_time
        efficiency = 100.0 * speedup / ideal_speedup
        p = ideal_speedup
        if math.isclose(p, 1.0):
            karp_flatt = 0.0
        else:
            karp_flatt = (1.0 / speedup - 1.0 / p) / (1.0 - 1.0 / p)

        row["seconds_per_iteration"] = median_time / iterations
        row["iterations_per_second"] = iterations / median_time
        row["relative_speedup"] = speedup
        row["ideal_relative_speedup"] = ideal_speedup
        row["relative_parallel_efficiency_pct"] = efficiency
        row["relative_karp_flatt"] = karp_flatt
        row["core_seconds_per_solve"] = ranks * median_time
        row["core_hours_per_solve"] = ranks * median_time / 3600.0
        setup = float(row["median_setup_time_s"])
        row["setup_fraction_pct"] = 100.0 * setup / (setup + median_time)

    output_dir = args.campaign_root / "analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = output_dir / "strong_scaling_summary.csv"

    fieldnames = [
        "ranks",
        "nodes",
        "job_id",
        "iterations",
        "measured_runs",
        "median_loop_time_s",
        "minimum_loop_time_s",
        "maximum_loop_time_s",
        "mean_loop_time_s",
        "stddev_loop_time_s",
        "median_setup_time_s",
        "median_elapsed_time_s",
        "seconds_per_iteration",
        "iterations_per_second",
        "relative_speedup",
        "ideal_relative_speedup",
        "relative_parallel_efficiency_pct",
        "relative_karp_flatt",
        "core_seconds_per_solve",
        "core_hours_per_solve",
        "setup_fraction_pct",
        "job_dir",
    ]

    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    markdown = output_dir / "strong_scaling_summary.md"
    with markdown.open("w", encoding="utf-8") as handle:
        handle.write("# CBS3D MPI strong-scaling summary\n\n")
        handle.write(f"Baseline: P{baseline_ranks}. Speedup and efficiency are relative to this measured baseline.\n\n")
        handle.write("| Ranks | Nodes | Median loop (s) | Speedup | Efficiency (%) | Iter/s | Core-hours | Rel. Karp-Flatt |\n")
        handle.write("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in rows:
            handle.write(
                f"| {int(row['ranks'])} | {int(row['nodes'])} | {row['median_loop_time_s']:.6f} "
                f"| {row['relative_speedup']:.4f} | {row['relative_parallel_efficiency_pct']:.2f} "
                f"| {row['iterations_per_second']:.3f} | {row['core_hours_per_solve']:.4f} "
                f"| {row['relative_karp_flatt']:.6f} |\n"
            )
        handle.write("\nThe relative Karp-Flatt value uses the smallest measured rank count as the reference. ")
        handle.write("A P40/P80 result is a two-point pilot; a defensible scaling curve requires the full P40--P1280 campaign.\n")

    plots_written = write_plots(rows, output_dir / "plots")

    print(f"Summary CSV : {summary_csv}")
    print(f"Summary MD  : {markdown}")
    if plots_written:
        print(f"Plots       : {output_dir / 'plots'}")
    else:
        print("Plots       : skipped (matplotlib not installed); run this script on a machine with matplotlib")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
