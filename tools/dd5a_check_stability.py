#!/usr/bin/env python3
"""Validate one CBS3D DD-5A long-run residual history and run log.

The checker is intentionally conservative and compatible with Sunbird's Python
3.6 interpreter. It verifies run completion, persistent PETSc setup, finite
iteration diagnostics, positive timesteps, bounded Krylov residuals and absence
of gross growth in continuity or velocity metrics.

This is a production-stability gate, not a proof of steady-state convergence.
"""

import argparse
import csv
import math
from pathlib import Path
import re
import statistics
import sys


REQUIRED_FLOAT_FIELDS = (
    "time",
    "dt",
    "u_rel",
    "u_norm",
    "u_abs",
    "v_rel",
    "v_norm",
    "v_abs",
    "w_rel",
    "w_norm",
    "w_abs",
    "p_rel",
    "p_norm",
    "p_abs",
    "T_rel",
    "T_norm",
    "T_abs",
    "velocity_rel_max",
    "continuity_rms",
    "continuity_max",
    "maximum_velocity",
    "maximum_velocity_correction",
    "cg_initial_l2",
    "cg_final_l2",
    "cg_relative_l2",
    "cg_max_abs",
    "iteration_wall_seconds",
)


def finite_float(text, row_number, field):
    try:
        value = float(text)
    except ValueError:
        raise ValueError(
            "row {} field {} is not numeric: {!r}".format(
                row_number,
                field,
                text,
            )
        )

    if not math.isfinite(value):
        raise ValueError(
            "row {} field {} is non-finite: {!r}".format(
                row_number,
                field,
                text,
            )
        )

    return value


def read_residuals(path, expected_iterations):
    if not path.is_file():
        raise FileNotFoundError("Residual CSV does not exist: {}".format(path))

    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)

    if reader.fieldnames is None:
        raise ValueError("Residual CSV has no header")

    required = {"iteration", "cg_iterations"}.union(REQUIRED_FLOAT_FIELDS)
    missing = sorted(required.difference(reader.fieldnames))
    if missing:
        raise ValueError(
            "Residual CSV is missing columns: {}".format(", ".join(missing))
        )

    if len(rows) != expected_iterations:
        raise ValueError(
            "Expected {} residual rows, found {}".format(
                expected_iterations,
                len(rows),
            )
        )

    parsed = []

    for row_number, row in enumerate(rows, start=1):
        try:
            iteration = int(row["iteration"])
            cg_iterations = int(row["cg_iterations"])
        except ValueError:
            raise ValueError(
                "row {} has invalid integer iteration data".format(row_number)
            )

        if iteration != row_number:
            raise ValueError(
                "Residual iteration sequence mismatch at row {}: {}".format(
                    row_number,
                    iteration,
                )
            )
        if cg_iterations < 0:
            raise ValueError(
                "row {} has negative CG iteration count".format(row_number)
            )

        values = {
            field: finite_float(row[field], row_number, field)
            for field in REQUIRED_FLOAT_FIELDS
        }
        values["iteration"] = iteration
        values["cg_iterations"] = cg_iterations
        parsed.append(values)

    return parsed


def extract_last_integer(text, label):
    pattern = re.compile(
        r"^\s*{}\s*:\s*(\d+)\s*$".format(re.escape(label)),
        re.MULTILINE,
    )
    matches = pattern.findall(text)
    return int(matches[-1]) if matches else None


def validate_log(path, expected_iterations, expected_ranks):
    if not path.is_file():
        raise FileNotFoundError("Run log does not exist: {}".format(path))

    text = path.read_text(encoding="utf-8", errors="replace")
    failures = []

    if re.search(r"\bERROR\b|failed to converge|MPI_ABORT", text, re.IGNORECASE):
        failures.append("run log contains ERROR, MPI_ABORT or convergence failure")

    if "CBS3D DD-4B/DD-4C PRODUCTION LOOP COMPLETE" not in text:
        failures.append("production-loop completion banner is missing")

    completed = extract_last_integer(text, "completed iterations")
    ranks = extract_last_integer(text, "MPI ranks")
    matrix_builds = extract_last_integer(text, "pressure matrix builds")
    amg_builds = extract_last_integer(text, "AMG hierarchy builds")

    if completed != expected_iterations:
        failures.append(
            "completed iterations is {}, expected {}".format(
                completed,
                expected_iterations,
            )
        )
    if ranks != expected_ranks:
        failures.append(
            "MPI ranks is {}, expected {}".format(ranks, expected_ranks)
        )
    if matrix_builds != 1:
        failures.append(
            "pressure matrix builds is {}, expected 1".format(matrix_builds)
        )
    if amg_builds != 1:
        failures.append(
            "AMG hierarchy builds is {}, expected 1".format(amg_builds)
        )

    return failures, {
        "completed_iterations": completed,
        "mpi_ranks": ranks,
        "pressure_matrix_builds": matrix_builds,
        "amg_hierarchy_builds": amg_builds,
    }


def percentile(values, fraction):
    ordered = sorted(values)
    if not ordered:
        return 0.0
    index = int(math.ceil(fraction * len(ordered))) - 1
    index = max(0, min(index, len(ordered) - 1))
    return ordered[index]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("residual_csv", type=Path)
    parser.add_argument("run_log", type=Path)
    parser.add_argument("--expected-iterations", type=int, default=100)
    parser.add_argument("--expected-ranks", type=int, default=40)
    parser.add_argument("--maximum-cg-relative", type=float, default=1.0e-3)
    args = parser.parse_args()

    if args.expected_iterations < 1:
        raise ValueError("expected-iterations must be positive")
    if args.expected_ranks < 2:
        raise ValueError("expected-ranks must be at least 2")
    if args.maximum_cg_relative <= 0.0:
        raise ValueError("maximum-cg-relative must be positive")

    rows = read_residuals(args.residual_csv, args.expected_iterations)
    log_failures, log_summary = validate_log(
        args.run_log,
        args.expected_iterations,
        args.expected_ranks,
    )

    failures = list(log_failures)

    dt_values = [row["dt"] for row in rows]
    cg_relative = [row["cg_relative_l2"] for row in rows]
    continuity = [row["continuity_rms"] for row in rows]
    continuity_max = [row["continuity_max"] for row in rows]
    velocity = [row["maximum_velocity"] for row in rows]
    correction = [row["maximum_velocity_correction"] for row in rows]
    wall = [row["iteration_wall_seconds"] for row in rows]
    cg_iterations = [row["cg_iterations"] for row in rows]

    if min(dt_values) <= 0.0:
        failures.append("one or more timesteps are non-positive")

    dt_ratio = max(dt_values) / max(min(dt_values), 1.0e-300)
    if dt_ratio > 1.01:
        failures.append(
            "timestep max/min ratio {:.6g} exceeds 1.01".format(dt_ratio)
        )

    maximum_cg_relative = max(cg_relative)
    if maximum_cg_relative > args.maximum_cg_relative:
        failures.append(
            "maximum CG relative residual {:.6e} exceeds {:.6e}".format(
                maximum_cg_relative,
                args.maximum_cg_relative,
            )
        )

    first_continuity = max(continuity[0], 1.0e-300)
    continuity_peak_ratio = max(continuity) / first_continuity
    continuity_final_ratio = continuity[-1] / first_continuity

    if continuity_peak_ratio > 100.0:
        failures.append(
            "continuity RMS peak/initial ratio {:.6g} exceeds 100".format(
                continuity_peak_ratio
            )
        )
    if continuity_final_ratio > 10.0:
        failures.append(
            "continuity RMS final/initial ratio {:.6g} exceeds 10".format(
                continuity_final_ratio
            )
        )

    first_velocity = max(velocity[0], 1.0e-300)
    velocity_peak_ratio = max(velocity) / first_velocity
    if velocity_peak_ratio > 10.0:
        failures.append(
            "maximum velocity peak/initial ratio {:.6g} exceeds 10".format(
                velocity_peak_ratio
            )
        )

    first_correction = max(correction[0], 1.0e-300)
    correction_peak_ratio = max(correction) / first_correction
    if correction_peak_ratio > 10.0:
        failures.append(
            "velocity-correction peak/initial ratio {:.6g} exceeds 10".format(
                correction_peak_ratio
            )
        )

    if min(wall) <= 0.0:
        failures.append("one or more iteration wall times are non-positive")

    print("DD-5A long-run stability report")
    print("  residual CSV              : {}".format(args.residual_csv))
    print("  run log                   : {}".format(args.run_log))
    print("  completed iterations      : {}".format(
        log_summary["completed_iterations"]
    ))
    print("  MPI ranks                 : {}".format(log_summary["mpi_ranks"]))
    print("  pressure matrix builds    : {}".format(
        log_summary["pressure_matrix_builds"]
    ))
    print("  AMG hierarchy builds      : {}".format(
        log_summary["amg_hierarchy_builds"]
    ))
    print("  dt min / max              : {:.8e} / {:.8e}".format(
        min(dt_values), max(dt_values)
    ))
    print("  CG iterations min / max   : {} / {}".format(
        min(cg_iterations), max(cg_iterations)
    ))
    print("  CG relative residual max  : {:.8e}".format(maximum_cg_relative))
    print("  continuity RMS first/final: {:.8e} / {:.8e}".format(
        continuity[0], continuity[-1]
    ))
    print("  continuity RMS peak       : {:.8e}".format(max(continuity)))
    print("  continuity max peak       : {:.8e}".format(max(continuity_max)))
    print("  Umax first/final          : {:.8e} / {:.8e}".format(
        velocity[0], velocity[-1]
    ))
    print("  Umax peak                 : {:.8e}".format(max(velocity)))
    print("  dUmax first/final         : {:.8e} / {:.8e}".format(
        correction[0], correction[-1]
    ))
    print("  iteration wall mean       : {:.8e} s".format(
        statistics.mean(wall)
    ))
    print("  iteration wall p95        : {:.8e} s".format(
        percentile(wall, 0.95)
    ))

    if failures:
        print()
        print("DD-5A long-run stability: FAIL")
        for failure in failures:
            print("  - {}".format(failure))
        return 1

    print()
    print("DD-5A long-run stability: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        raise SystemExit(2)
