#!/usr/bin/env python3
"""
Automatically refresh CBS3D++_SI flat-plate SA post-processing outputs.

This watcher is intentionally a thin wrapper around postprocess_flatplate_sa.py.
It does not duplicate the numerical extraction logic.  It simply watches a case
output directory and reruns the post-processor whenever a new VTU file or
residual CSV update is detected.

Typical use from the repository root:

    python tools/watch_flatplate_sa_outputs.py \
        output/input/flat_plate_sa_3d_corrected/flat_plate_sa_3d_corrected \
        --rho 1.0 --mu 1.0e-6 --u-inf 1.0 \
        --stations 0.97 1.90 \
        --interval 60

The generated CSV and PNG files are written by postprocess_flatplate_sa.py into:

    <case-output-dir>/post_flatplate_sa
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable, List, Optional, Tuple


ROOT = Path(__file__).resolve().parents[1]
POSTPROCESS = ROOT / "tools" / "postprocess_flatplate_sa.py"


def latest_mtime(paths: Iterable[Path]) -> float:
    value = 0.0
    for path in paths:
        try:
            value = max(value, path.stat().st_mtime)
        except FileNotFoundError:
            continue
    return value


def watched_files(output_dir: Path) -> List[Path]:
    files: List[Path] = []
    files.extend(output_dir.rglob("*.vtu"))
    files.extend(output_dir.rglob("*_residuals.csv"))
    return files


def latest_vtu(output_dir: Path) -> Optional[Path]:
    files = list(output_dir.rglob("*.vtu"))
    if not files:
        return None
    return max(files, key=lambda p: p.stat().st_mtime)


def run_postprocessor(args: argparse.Namespace) -> int:
    cmd = [
        sys.executable,
        str(POSTPROCESS),
        str(args.output_dir),
        "--rho",
        str(args.rho),
        "--mu",
        str(args.mu),
        "--u-inf",
        str(args.u_inf),
        "--wall-bc",
        str(args.wall_bc),
        "--station-dx",
        str(args.station_dx),
        "--cf-bins",
        str(args.cf_bins),
    ]

    if args.stations:
        cmd.append("--stations")
        cmd.extend(str(x) for x in args.stations)

    if args.no_plots:
        cmd.append("--no-plots")

    print("[post] " + " ".join(cmd), flush=True)
    completed = subprocess.run(cmd, cwd=ROOT)
    return completed.returncode


def run_once_when_ready(args: argparse.Namespace) -> int:
    if not args.output_dir.exists():
        raise SystemExit(f"Output directory does not exist: {args.output_dir}")

    vtu = latest_vtu(args.output_dir)
    if vtu is None:
        print(f"[wait] No VTU files found yet under {args.output_dir}", flush=True)
        return 1

    print(f"[data] Latest VTU: {vtu}", flush=True)
    return run_postprocessor(args)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Watch a CBS3D++_SI flat-plate SA output directory and refresh plots automatically."
    )
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--rho", type=float, required=True)
    parser.add_argument("--mu", type=float, required=True)
    parser.add_argument("--u-inf", dest="u_inf", type=float, required=True)
    parser.add_argument("--wall-bc", type=int, default=530)
    parser.add_argument("--stations", type=float, nargs="*", default=[0.97, 1.90])
    parser.add_argument("--station-dx", type=float, default=0.02)
    parser.add_argument("--cf-bins", type=int, default=100)
    parser.add_argument("--interval", type=float, default=60.0)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    args.output_dir = args.output_dir.resolve()

    if args.once:
        raise SystemExit(run_once_when_ready(args))

    print("[watch] CBS3D++_SI flat-plate SA post-processing watcher", flush=True)
    print(f"[watch] Output directory : {args.output_dir}", flush=True)
    print(f"[watch] Refresh interval : {args.interval:g} s", flush=True)
    print("[watch] Press Ctrl+C to stop.\n", flush=True)

    last_seen: Optional[Tuple[float, int]] = None

    while True:
        files = watched_files(args.output_dir)
        stamp = latest_mtime(files)
        count = len(files)
        current = (stamp, count)

        if count == 0:
            print(f"[wait] No VTU/residual files found yet under {args.output_dir}", flush=True)
        elif current != last_seen:
            vtu = latest_vtu(args.output_dir)
            if vtu is not None:
                print(f"[data] Detected output update. Latest VTU: {vtu.name}", flush=True)
                code = run_postprocessor(args)
                if code != 0:
                    print(
                        "[warn] Post-processing failed. This can happen if the solver is still writing a VTU file. "
                        "The watcher will retry at the next interval.",
                        flush=True,
                    )
                else:
                    print("[ok] Post-processing refreshed.\n", flush=True)
            last_seen = current
        else:
            print("[idle] No new solver output detected.", flush=True)

        time.sleep(args.interval)


if __name__ == "__main__":
    main()
