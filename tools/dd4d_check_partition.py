#!/usr/bin/env python3
"""Validate a CBS3D MPI partition tree before a rank-count run.

The checker verifies that every requested rank has the complete local case and
that each .mpi metadata header reports the expected rank and communicator size.
It deliberately performs no mesh parsing and is therefore fast enough to run on
a login node.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


REQUIRED_EXTENSIONS = ("plt", "bco", "par", "mpi", "material", "matprop")


def rank_text(rank: int) -> str:
    return f"{rank:04d}"


def parse_partition_header(path: Path) -> tuple[int, int, int, int, int, int]:
    """Return partition_id, mpi_rank, mpi_size, global nelem/npoin/nboun."""
    with path.open("r", encoding="utf-8") as stream:
        first = stream.readline().split()
        second = stream.readline().split()
        third = stream.readline().split()

    if len(first) != 2 or first[0] != "CBS3D_MPI_PARTITION":
        raise ValueError(f"{path}: invalid CBS3D_MPI_PARTITION header")

    if int(first[1]) != 1:
        raise ValueError(f"{path}: unsupported metadata version {first[1]}")

    if len(second) != 4 or second[0] != "PARTITION":
        raise ValueError(f"{path}: invalid PARTITION header")

    if len(third) != 4 or third[0] != "GLOBAL":
        raise ValueError(f"{path}: invalid GLOBAL header")

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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("partition_root", type=Path)
    parser.add_argument("case_name")
    parser.add_argument("expected_ranks", type=int)
    args = parser.parse_args()

    root: Path = args.partition_root.resolve()
    case_name: str = args.case_name
    expected_ranks: int = args.expected_ranks

    if expected_ranks < 2:
        raise ValueError("expected_ranks must be at least 2")

    if not root.is_dir():
        raise FileNotFoundError(f"Partition root does not exist: {root}")

    reference_global: tuple[int, int, int] | None = None
    total_missing = 0
    total_header_errors = 0

    print(f"Partition root : {root}")
    print(f"Case           : {case_name}")
    print(f"Expected ranks : {expected_ranks}")
    print()

    for rank in range(expected_ranks):
        text = rank_text(rank)
        rank_dir = root / f"rank_{text}"
        base = rank_dir / f"{case_name}_rank_{text}"

        missing = [
            str(base.with_suffix(f".{extension}"))
            for extension in REQUIRED_EXTENSIONS
            if not base.with_suffix(f".{extension}").is_file()
        ]

        if missing:
            total_missing += len(missing)
            print(f"rank {text}: FAIL")
            for path in missing:
                print(f"  missing: {path}")
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
                    f"partition_id={partition_id}, expected {rank}"
                )

            if metadata_rank != rank:
                raise ValueError(
                    f"mpi_rank={metadata_rank}, expected {rank}"
                )

            if metadata_size != expected_ranks:
                raise ValueError(
                    f"mpi_size={metadata_size}, expected {expected_ranks}"
                )

            global_counts = (global_nelem, global_npoin, global_nboun)

            if reference_global is None:
                reference_global = global_counts
            elif global_counts != reference_global:
                raise ValueError(
                    f"GLOBAL={global_counts}, expected {reference_global}"
                )

        except (OSError, ValueError) as error:
            total_header_errors += 1
            print(f"rank {text}: FAIL")
            print(f"  metadata: {error}")
            continue

        print(
            f"rank {text}: PASS  "
            f"GLOBAL={global_nelem}/{global_npoin}/{global_nboun}"
        )

    print()

    if total_missing or total_header_errors:
        print(
            "DD-4D partition validation: FAIL  "
            f"missing_files={total_missing} "
            f"metadata_errors={total_header_errors}"
        )
        return 1

    assert reference_global is not None
    print(
        "DD-4D partition validation: PASS  "
        f"ranks={expected_ranks} "
        f"global_nelem={reference_global[0]} "
        f"global_npoin={reference_global[1]} "
        f"global_nboun={reference_global[2]}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:  # noqa: BLE001 - command-line diagnostic
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
