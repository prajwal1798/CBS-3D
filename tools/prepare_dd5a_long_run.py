#!/usr/bin/env python3
"""Prepare an isolated CBS3D DD-5A long-run partition tree.

The script is compatible with the Python 3.6 interpreter on Sunbird. It creates
an efficient hard-linked copy of an existing partition tree, then atomically
replaces only the rank-local .par files. The source partition is never edited.

DD-5A policy:

* run the requested number of iterations;
* retain one residual row per iteration;
* print console diagnostics at a controlled cadence;
* disable live plotting and VTU/PVTU/PVD output;
* prevent steady-state early termination before the requested final iteration.

Disabling field output here is deliberate. DD-5A isolates numerical and MPI
long-run stability from the separate DD-5B scalable-output implementation.
"""

import argparse
import hashlib
import os
from pathlib import Path
import shutil
import sys


REQUIRED_EXTENSIONS = ("plt", "bco", "par", "mpi", "material", "matprop")
TIMESTEP_LABEL = "ntime, transient_on, dtfixed, dtfix, iwrite"
OUTPUT_LABEL = "cbs3d output/monitor controls:"


def rank_text(rank):
    return "{:04d}".format(rank)


def link_or_copy(source, destination):
    """Hard-link files when possible and copy only when linking is unavailable."""
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)


def next_data_line(lines, label_index, path):
    for index in range(label_index + 1, len(lines)):
        stripped = lines[index].strip()
        if not stripped:
            continue
        if stripped.startswith("#") or stripped.startswith("!"):
            continue
        return index

    raise ValueError(
        "{}: no data line follows label at line {}".format(
            path,
            label_index + 1,
        )
    )


def find_unique_label(lines, needle, path):
    matches = []
    lowered = needle.lower()

    for index, line in enumerate(lines):
        if lowered in line.lower():
            matches.append(index)

    if len(matches) != 1:
        raise ValueError(
            "{}: expected one label containing {!r}, found {}".format(
                path,
                needle,
                len(matches),
            )
        )

    return matches[0]


def replace_tokens_preserving_indent(original_line, tokens):
    indent_length = len(original_line) - len(original_line.lstrip(" \t"))
    indent = original_line[:indent_length]
    newline = "\n" if original_line.endswith("\n") else ""
    return indent + " ".join(tokens) + newline


def update_parameter_file(path, iterations, console_every):
    with path.open("r", encoding="utf-8") as stream:
        lines = stream.readlines()

    timestep_label_index = find_unique_label(lines, TIMESTEP_LABEL, path)
    timestep_data_index = next_data_line(lines, timestep_label_index, path)
    timestep_tokens = lines[timestep_data_index].split()

    if len(timestep_tokens) < 5:
        raise ValueError(
            "{}: timestep control requires at least five values".format(path)
        )

    timestep_tokens[0] = str(iterations)
    lines[timestep_data_index] = replace_tokens_preserving_indent(
        lines[timestep_data_index],
        timestep_tokens,
    )

    output_label_index = find_unique_label(lines, OUTPUT_LABEL, path)
    output_data_index = next_data_line(lines, output_label_index, path)
    output_tokens = lines[output_data_index].split()

    if len(output_tokens) < 9:
        raise ValueError(
            "{}: output/monitor control requires at least nine values".format(
                path
            )
        )

    # residual_log_enabled residual_log_every console_log_every
    # live_residual_plot vtu_output_enabled vtu_output_every_iterations
    # vtu_output_every_sim_time write_boundary_debug_arrays
    # steady_min_iterations [future optional values...]
    output_tokens[0] = "1"
    output_tokens[1] = "1"
    output_tokens[2] = str(console_every)
    output_tokens[3] = "0"
    output_tokens[4] = "0"
    output_tokens[5] = str(iterations)
    output_tokens[6] = "0"
    output_tokens[8] = str(iterations + 1)

    lines[output_data_index] = replace_tokens_preserving_indent(
        lines[output_data_index],
        output_tokens,
    )

    temporary = path.with_name(path.name + ".dd5a.tmp")

    with temporary.open("w", encoding="utf-8") as stream:
        stream.writelines(lines)
        stream.flush()
        os.fsync(stream.fileno())

    os.replace(str(temporary), str(path))

    return {
        "timestep_line": timestep_data_index + 1,
        "output_line": output_data_index + 1,
        "timestep_values": " ".join(timestep_tokens),
        "output_values": " ".join(output_tokens),
    }


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            block = stream.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def validate_rank_files(root, case_name, ranks):
    missing = []

    for rank in range(ranks):
        text = rank_text(rank)
        base = root / "rank_{}".format(text) / "{}_rank_{}".format(
            case_name,
            text,
        )

        for extension in REQUIRED_EXTENSIONS:
            path = base.with_suffix(".{}".format(extension))
            if not path.is_file():
                missing.append(str(path))

    if missing:
        raise FileNotFoundError(
            "Prepared tree is incomplete; missing {} files:\n{}".format(
                len(missing),
                "\n".join(missing[:40]),
            )
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_root", type=Path)
    parser.add_argument("destination_root", type=Path)
    parser.add_argument("case_name")
    parser.add_argument("ranks", type=int)
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--console-every", type=int, default=10)
    args = parser.parse_args()

    source = args.source_root.resolve()
    destination = args.destination_root.resolve()

    if not source.is_dir():
        raise FileNotFoundError(
            "Source partition root does not exist: {}".format(source)
        )
    if destination.exists():
        raise FileExistsError(
            "Destination already exists: {}".format(destination)
        )
    if args.ranks < 2:
        raise ValueError("ranks must be at least 2")
    if args.iterations < 1:
        raise ValueError("iterations must be at least 1")
    if args.console_every < 1:
        raise ValueError("console-every must be at least 1")

    validate_rank_files(source, args.case_name, args.ranks)

    print("Creating isolated DD-5A partition tree")
    print("  source      : {}".format(source))
    print("  destination : {}".format(destination))
    print("  ranks       : {}".format(args.ranks))
    print("  iterations  : {}".format(args.iterations))
    print("  VTU output  : disabled")

    shutil.copytree(
        str(source),
        str(destination),
        copy_function=link_or_copy,
        symlinks=False,
    )

    parameter_hashes = set()
    first_summary = None

    for rank in range(args.ranks):
        text = rank_text(rank)
        parameter_path = (
            destination
            / "rank_{}".format(text)
            / "{}_rank_{}.par".format(args.case_name, text)
        )
        summary = update_parameter_file(
            parameter_path,
            args.iterations,
            args.console_every,
        )
        parameter_hashes.add(sha256(parameter_path))
        if first_summary is None:
            first_summary = summary

    validate_rank_files(destination, args.case_name, args.ranks)

    if len(parameter_hashes) != 1:
        raise RuntimeError(
            "Prepared rank-local parameter files are not identical"
        )

    print()
    print("DD-5A partition preparation: PASS")
    print("  parameter SHA-256 : {}".format(next(iter(parameter_hashes))))
    print("  timestep line     : {}".format(first_summary["timestep_line"]))
    print("  timestep values   : {}".format(first_summary["timestep_values"]))
    print("  output line       : {}".format(first_summary["output_line"]))
    print("  output values     : {}".format(first_summary["output_values"]))
    print("  source unchanged  : YES")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        raise SystemExit(2)
