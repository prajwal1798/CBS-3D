#!/usr/bin/env python3
"""Validate a CBS3D MPI partition tree before a rank-count run.

This script is intentionally compatible with the system Python 3.6 available on
Sunbird. It checks that every requested rank has the complete local case and
that each .mpi metadata header reports the expected rank, communicator size and
common global mesh counts.
"""

import argparse
from pathlib import Path
import sys


REQUIRED_EXTENSIONS = ("plt", "bco", "par", "mpi", "material", "matprop")


def rank_text(rank):
    return "{:04d}".format(rank)


def parse_partition_header(path):
    """Return partition_id, mpi_rank, mpi_size, global nelem/npoin/nboun."""
    with path.open("r", encoding="utf-8") as stream:
        first = stream.readline().split()
        second = stream.readline().split()
        third = stream.readline().split()

    if len(first) != 2 or first[0] != "CBS3D_MPI_PARTITION":
        raise ValueError("{}: invalid CBS3D_MPI_PARTITION header".format(path))

    if int(first[1]) != 1:
        raise ValueError(
            "{}: unsupported metadata version {}".format(path, first[1])
        )

    if len(second) != 4 or second[0] != "PARTITION":
        raise ValueError("{}: invalid PARTITION header".format(path))

    if len(third) != 4 or third[0] != "GLOBAL":
        raise ValueError("{}: invalid GLOBAL header".format(path))

    partition_id, mpi_rank, mpi_size = map(int, second[1:])
    global_nelem, global_npoin, global_nboun = map(int, third[1:])

    return (
        partition_id,
        mpi_rank,
        mpi_size,
        global_nelem,
        global_npoin,
        global_nboun,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("partition_root", type=Path)
    parser.add_argument("case_name")
    parser.add_argument("expected_ranks", type=int)
    args = parser.parse_args()

    root = args.partition_root.resolve()
    case_name = args.case_name
    expected_ranks = args.expected_ranks

    if expected_ranks < 2:
        raise ValueError("expected_ranks must be at least 2")

    if not root.is_dir():
        raise FileNotFoundError(
            "Partition root does not exist: {}".format(root)
        )

    reference_global = None
    total_missing = 0
    total_header_errors = 0

    print("Partition root : {}".format(root))
    print("Case           : {}".format(case_name))
    print("Expected ranks : {}".format(expected_ranks))
    print()

    for rank in range(expected_ranks):
        text = rank_text(rank)
        rank_dir = root / "rank_{}".format(text)
        base = rank_dir / "{}_rank_{}".format(case_name, text)

        missing = []
        for extension in REQUIRED_EXTENSIONS:
            candidate = base.with_suffix(".{}".format(extension))
            if not candidate.is_file():
                missing.append(str(candidate))

        if missing:
            total_missing += len(missing)
            print("rank {}: FAIL".format(text))
            for path in missing:
                print("  missing: {}".format(path))
            continue

        metadata_path = base.with_suffix(".mpi")

        try:
            (
                partition_id,
                metadata_rank,
                metadata_size,
                global_nelem,
                global_npoin,
                global_nboun,
            ) = parse_partition_header(metadata_path)

            if partition_id != rank:
                raise ValueError(
                    "partition_id={}, expected {}".format(partition_id, rank)
                )

            if metadata_rank != rank:
                raise ValueError(
                    "mpi_rank={}, expected {}".format(metadata_rank, rank)
                )

            if metadata_size != expected_ranks:
                raise ValueError(
                    "mpi_size={}, expected {}".format(
                        metadata_size,
                        expected_ranks,
                    )
                )

            global_counts = (global_nelem, global_npoin, global_nboun)

            if reference_global is None:
                reference_global = global_counts
            elif global_counts != reference_global:
                raise ValueError(
                    "GLOBAL={}, expected {}".format(
                        global_counts,
                        reference_global,
                    )
                )

        except (OSError, ValueError) as error:
            total_header_errors += 1
            print("rank {}: FAIL".format(text))
            print("  metadata: {}".format(error))
            continue

        print(
            "rank {}: PASS  GLOBAL={}/{}/{}".format(
                text,
                global_nelem,
                global_npoin,
                global_nboun,
            )
        )

    print()

    if total_missing or total_header_errors:
        print(
            "DD-4D partition validation: FAIL  "
            "missing_files={} metadata_errors={}".format(
                total_missing,
                total_header_errors,
            )
        )
        return 1

    if reference_global is None:
        raise RuntimeError("No valid partition metadata was found")

    print(
        "DD-4D partition validation: PASS  "
        "ranks={} global_nelem={} global_npoin={} global_nboun={}".format(
            expected_ranks,
            reference_global[0],
            reference_global[1],
            reference_global[2],
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        raise SystemExit(2)
