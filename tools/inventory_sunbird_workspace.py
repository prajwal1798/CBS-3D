#!/usr/bin/env python3
"""Inventory a CBS3D Sunbird workspace without modifying or deleting files.

The report is intended to precede consolidation and cleanup. It records the
size, likely role, Git state and selected file counts for every top-level item
under a requested scratch directory.

The script deliberately uses subprocess ``cwd=`` rather than ``git -C`` because
older Git installations on Sunbird do not support the ``-C`` option.
"""

from __future__ import print_function

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime


INPUT_SUFFIXES = {
    ".bco",
    ".material",
    ".matprop",
    ".mpi",
    ".par",
    ".plt",
}

RESULT_SUFFIXES = {
    ".csv",
    ".err",
    ".log",
    ".out",
    ".pvd",
    ".pvtu",
    ".vtu",
}


def run(command, cwd=None):
    try:
        process = subprocess.Popen(
            command,
            cwd=str(cwd) if cwd is not None else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = process.communicate()
        return process.returncode, stdout.strip(), stderr.strip()
    except OSError as error:
        return 127, "", str(error)


def directory_size_bytes(path):
    code, output, _ = run(["du", "-sb", str(path)])

    if code == 0 and output:
        try:
            return int(output.split()[0])
        except (ValueError, IndexError):
            pass

    code, output, _ = run(["du", "-sk", str(path)])

    if code == 0 and output:
        try:
            return int(output.split()[0]) * 1024
        except (ValueError, IndexError):
            pass

    if path.is_file():
        try:
            return path.stat().st_size
        except OSError:
            return 0

    return 0


def human_size(size):
    value = float(size)

    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return "{:.2f} {}".format(value, unit)
        value /= 1024.0

    return "{:.2f} TiB".format(value)


def git_metadata(path):
    metadata = {
        "is_git_repository": False,
        "git_branch": "",
        "git_commit": "",
        "git_remote": "",
        "git_dirty": False,
        "git_status_lines": 0,
        "git_error": "",
    }

    code, _, error = run(
        ["git", "rev-parse", "--is-inside-work-tree"],
        cwd=path,
    )

    if code != 0:
        metadata["git_error"] = error
        return metadata

    metadata["is_git_repository"] = True

    code, branch, _ = run(
        ["git", "symbolic-ref", "--short", "-q", "HEAD"],
        cwd=path,
    )
    metadata["git_branch"] = branch if code == 0 else "DETACHED"

    code, commit, error = run(
        ["git", "rev-parse", "HEAD"],
        cwd=path,
    )
    if code == 0:
        metadata["git_commit"] = commit
    else:
        metadata["git_error"] = error

    code, remote, _ = run(
        ["git", "config", "--get", "remote.origin.url"],
        cwd=path,
    )
    if code == 0:
        metadata["git_remote"] = remote

    code, status, error = run(
        ["git", "status", "--porcelain"],
        cwd=path,
    )
    if code == 0:
        lines = [line for line in status.splitlines() if line.strip()]
        metadata["git_dirty"] = bool(lines)
        metadata["git_status_lines"] = len(lines)
    elif not metadata["git_error"]:
        metadata["git_error"] = error

    return metadata


def inspect_tree(path, maximum_files):
    counts = {
        "input_files": 0,
        "result_files": 0,
        "slurm_files": 0,
        "cmake_cache_files": 0,
        "total_files_sampled": 0,
        "scan_truncated": False,
    }

    if path.is_file():
        suffix = path.suffix.lower()
        counts["input_files"] = int(suffix in INPUT_SUFFIXES)
        counts["result_files"] = int(suffix in RESULT_SUFFIXES)
        counts["slurm_files"] = int(suffix == ".slurm")
        counts["cmake_cache_files"] = int(path.name == "CMakeCache.txt")
        counts["total_files_sampled"] = 1
        return counts

    try:
        for root, directories, files in os.walk(str(path)):
            directories[:] = [
                name
                for name in directories
                if name not in {".git", "__pycache__"}
            ]

            for name in files:
                counts["total_files_sampled"] += 1
                suffix = Path(name).suffix.lower()

                if suffix in INPUT_SUFFIXES:
                    counts["input_files"] += 1
                if suffix in RESULT_SUFFIXES:
                    counts["result_files"] += 1
                if suffix == ".slurm":
                    counts["slurm_files"] += 1
                if name == "CMakeCache.txt":
                    counts["cmake_cache_files"] += 1

                if counts["total_files_sampled"] >= maximum_files:
                    counts["scan_truncated"] = True
                    return counts
    except OSError:
        counts["scan_truncated"] = True

    return counts


def classify(path, git, counts):
    name = path.name.lower()

    if path.is_file():
        return "file"

    if git["is_git_repository"]:
        return "source-repository"

    if counts["cmake_cache_files"] > 0 or name.startswith("build"):
        return "build-directory"

    if "backup" in name or "recovery" in name or "archive" in name:
        return "backup-or-recovery"

    if "run" in name or counts["result_files"] > 20:
        return "run-or-results"

    if "benchmark" in name:
        return "benchmark-workspace"

    if "partition" in name:
        return "partition-data"

    if counts["input_files"] > 0:
        return "case-or-input-data"

    if name in {"projects", "project"}:
        return "project-container"

    return "unclassified-directory"


def cleanup_candidate(classification, git, counts):
    if classification == "build-directory":
        return "REVIEW_FOR_DELETION_REPRODUCIBLE_BUILD"

    if classification == "source-repository" and git["git_dirty"]:
        return "KEEP_UNCOMMITTED_SOURCE"

    if classification == "source-repository":
        return "COMPARE_COMMIT_BEFORE_ARCHIVE"

    if classification in {
        "run-or-results",
        "benchmark-workspace",
        "partition-data",
        "case-or-input-data",
        "backup-or-recovery",
    }:
        return "KEEP_UNTIL_AUTHORITATIVE_COPY_IDENTIFIED"

    if counts["result_files"] > 0 or counts["input_files"] > 0:
        return "KEEP_PENDING_CONTENT_REVIEW"

    return "REVIEW_MANUALLY"


def inspect_item(path, maximum_files):
    size = directory_size_bytes(path)
    git = git_metadata(path) if path.is_dir() else git_metadata(path.parent)

    if path.is_file():
        git = {
            "is_git_repository": False,
            "git_branch": "",
            "git_commit": "",
            "git_remote": "",
            "git_dirty": False,
            "git_status_lines": 0,
            "git_error": "",
        }

    counts = inspect_tree(path, maximum_files)
    classification = classify(path, git, counts)

    record = {
        "name": path.name,
        "path": str(path),
        "item_type": "directory" if path.is_dir() else "file",
        "classification": classification,
        "size_bytes": size,
        "size_human": human_size(size),
        "modified_time": "",
        "cleanup_status": cleanup_candidate(classification, git, counts),
    }

    try:
        record["modified_time"] = datetime.fromtimestamp(
            path.stat().st_mtime
        ).isoformat()
    except OSError:
        pass

    record.update(git)
    record.update(counts)
    return record


def write_csv(path, records):
    fields = [
        "name",
        "path",
        "item_type",
        "classification",
        "size_bytes",
        "size_human",
        "modified_time",
        "cleanup_status",
        "is_git_repository",
        "git_branch",
        "git_commit",
        "git_remote",
        "git_dirty",
        "git_status_lines",
        "git_error",
        "input_files",
        "result_files",
        "slurm_files",
        "cmake_cache_files",
        "total_files_sampled",
        "scan_truncated",
    ]

    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Inventory CBS3D scratch directories without deleting files"
    )
    parser.add_argument(
        "--root",
        required=True,
        help="Top-level Sunbird workspace, for example /scratch/s.2337862",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Report directory; defaults to <root>/CBS3D_inventory_<timestamp>",
    )
    parser.add_argument(
        "--maximum-files-per-item",
        type=int,
        default=200000,
        help="Safety limit for recursive file counting",
    )
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    root = Path(arguments.root).expanduser().resolve()

    if not root.is_dir():
        raise SystemExit("Workspace root does not exist: {}".format(root))

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    if arguments.output_dir is None:
        output_dir = root / ("CBS3D_inventory_" + timestamp)
    else:
        output_dir = Path(arguments.output_dir).expanduser().resolve()

    output_dir.mkdir(parents=True, exist_ok=False)

    records = []

    for path in sorted(root.iterdir(), key=lambda item: item.name.lower()):
        if path == output_dir:
            continue

        print("Inspecting {}".format(path))
        records.append(
            inspect_item(path, arguments.maximum_files_per_item)
        )

    csv_path = output_dir / "workspace_inventory.csv"
    json_path = output_dir / "workspace_inventory.json"

    write_csv(csv_path, records)

    with json_path.open("w", encoding="utf-8") as stream:
        json.dump(
            {
                "generated_at": datetime.now().isoformat(),
                "root": str(root),
                "records": records,
            },
            stream,
            indent=2,
            sort_keys=True,
        )
        stream.write("\n")

    print()
    print("Inventory complete")
    print("CSV : {}".format(csv_path))
    print("JSON: {}".format(json_path))
    print()
    print("No files or directories were deleted.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
