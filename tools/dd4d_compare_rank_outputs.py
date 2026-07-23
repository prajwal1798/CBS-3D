#!/usr/bin/env python3
"""Compare CBS3D distributed solutions produced with two rank counts.

The comparison is field-level rather than image-level. It reads the final PVTU
state from each run, retains only owner copies, matches nodes by global_node_id,
and reports absolute and relative L2/max differences for u, v, w, pressure and
temperature. The script also compares the rank-zero residual CSV histories.

The implementation is compatible with the system Python 3.6 on Sunbird. VTU
pieces are scanned line by line and are never loaded as XML trees.
"""

import argparse
import csv
import json
import math
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET


STEP_PATTERN = re.compile(r"_step_(\d+)\.pvtu$")
FIELD_NAMES = ("u", "v", "w", "pressure", "temperature")
CSV_FIELDS = (
    "dt",
    "u_rel",
    "v_rel",
    "w_rel",
    "p_rel",
    "T_rel",
    "continuity_rms",
    "continuity_max",
    "maximum_velocity",
    "maximum_velocity_correction",
    "cg_initial_l2",
    "cg_final_l2",
    "cg_relative_l2",
)


class FieldTolerance(object):
    __slots__ = ("absolute", "relative")

    def __init__(self, absolute, relative):
        self.absolute = float(absolute)
        self.relative = float(relative)


DEFAULT_TOLERANCES = {
    "u": FieldTolerance(absolute=1.0e-9, relative=1.0e-5),
    "v": FieldTolerance(absolute=1.0e-9, relative=1.0e-5),
    "w": FieldTolerance(absolute=1.0e-9, relative=1.0e-5),
    "pressure": FieldTolerance(absolute=1.0e-5, relative=1.0e-5),
    "temperature": FieldTolerance(absolute=1.0e-7, relative=1.0e-8),
}


def final_pvtu(output_dir, case_name):
    candidates = []

    for path in output_dir.glob("{}_step_*.pvtu".format(case_name)):
        match = STEP_PATTERN.search(path.name)
        if match:
            candidates.append((int(match.group(1)), path))

    if not candidates:
        raise FileNotFoundError(
            "No {}_step_*.pvtu files in {}".format(case_name, output_dir)
        )

    return max(candidates, key=lambda item: item[0])


def pvtu_piece_paths(pvtu_path):
    tree = ET.parse(str(pvtu_path))
    pieces = []

    for element in tree.findall(".//Piece"):
        source = element.attrib.get("Source")
        if not source:
            raise ValueError(
                "{}: Piece has no Source attribute".format(pvtu_path)
            )

        path = pvtu_path.parent / source
        if not path.is_file():
            raise FileNotFoundError(
                "{}: referenced VTU piece is missing: {}".format(
                    pvtu_path,
                    path,
                )
            )
        pieces.append(path)

    if not pieces:
        raise ValueError("{}: contains no Piece elements".format(pvtu_path))

    return pieces


def read_data_array(path, name, converter):
    marker = 'Name="{}"'.format(name)
    values = []
    in_array = False

    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not in_array:
                if "<DataArray" in line and marker in line:
                    in_array = True
                continue

            if "</DataArray>" in line:
                return values

            for token in line.split():
                try:
                    values.append(converter(token))
                except ValueError as error:
                    raise ValueError(
                        "{}:{}: invalid value in {}: {}".format(
                            path,
                            line_number,
                            name,
                            token,
                        )
                    ) from error

    raise ValueError(
        "{}: DataArray {!r} was not terminated".format(path, name)
    )


def load_owned_solution(output_dir, case_name):
    iteration, pvtu = final_pvtu(output_dir, case_name)
    pieces = pvtu_piece_paths(pvtu)

    # Mutable record order: u, v, w, pressure, temperature.
    records = {}

    for piece_index, piece in enumerate(pieces, start=1):
        print(
            "  reading piece {}/{}: {}".format(
                piece_index,
                len(pieces),
                piece.name,
            ),
            flush=True,
        )

        global_ids = read_data_array(piece, "global_node_id", int)
        is_owned = read_data_array(piece, "is_owned", int)

        if len(global_ids) != len(is_owned):
            raise ValueError(
                "{}: global_node_id/is_owned length mismatch".format(piece)
            )

        owned_indices = [
            index for index, flag in enumerate(is_owned) if flag != 0
        ]

        for index in owned_indices:
            global_id = global_ids[index]
            if global_id in records:
                raise ValueError(
                    "Duplicate owned global node {} in {}".format(
                        global_id,
                        piece,
                    )
                )
            records[global_id] = [math.nan] * 5

        pressure = read_data_array(piece, "pressure", float)
        temperature = read_data_array(piece, "temperature", float)
        velocity = read_data_array(piece, "velocity", float)

        node_count = len(global_ids)

        if len(pressure) != node_count:
            raise ValueError("{}: pressure length mismatch".format(piece))
        if len(temperature) != node_count:
            raise ValueError("{}: temperature length mismatch".format(piece))
        if len(velocity) != 3 * node_count:
            raise ValueError("{}: velocity length mismatch".format(piece))

        for index in owned_indices:
            global_id = global_ids[index]
            offset = 3 * index
            record = records[global_id]
            record[0] = velocity[offset]
            record[1] = velocity[offset + 1]
            record[2] = velocity[offset + 2]
            record[3] = pressure[index]
            record[4] = temperature[index]

    immutable_records = {}

    for global_id, record in records.items():
        if not all(math.isfinite(value) for value in record):
            raise ValueError(
                "Non-finite or incomplete field record at global node {}".format(
                    global_id
                )
            )
        immutable_records[global_id] = tuple(record)

    return iteration, len(pieces), immutable_records


def compare_fields(reference, candidate, tolerances):
    reference_ids = set(reference)
    candidate_ids = set(candidate)

    missing = sorted(reference_ids - candidate_ids)
    extra = sorted(candidate_ids - reference_ids)

    if missing or extra:
        raise ValueError(
            "Owned global-node sets differ: missing={}, extra={}".format(
                len(missing),
                len(extra),
            )
        )

    accumulators = {}
    for name in FIELD_NAMES:
        accumulators[name] = {
            "difference_squared": 0.0,
            "reference_squared": 0.0,
            "candidate_squared": 0.0,
            "maximum_absolute_difference": 0.0,
            "maximum_reference_absolute": 0.0,
        }

    for global_id in reference_ids:
        reference_record = reference[global_id]
        candidate_record = candidate[global_id]

        for field_index, name in enumerate(FIELD_NAMES):
            reference_value = reference_record[field_index]
            candidate_value = candidate_record[field_index]
            difference = candidate_value - reference_value
            accumulator = accumulators[name]

            accumulator["difference_squared"] += difference * difference
            accumulator["reference_squared"] += (
                reference_value * reference_value
            )
            accumulator["candidate_squared"] += (
                candidate_value * candidate_value
            )
            accumulator["maximum_absolute_difference"] = max(
                accumulator["maximum_absolute_difference"],
                abs(difference),
            )
            accumulator["maximum_reference_absolute"] = max(
                accumulator["maximum_reference_absolute"],
                abs(reference_value),
            )

    all_passed = True
    report = {}

    for name in FIELD_NAMES:
        accumulator = accumulators[name]
        difference_l2 = math.sqrt(accumulator["difference_squared"])
        reference_l2 = math.sqrt(accumulator["reference_squared"])
        candidate_l2 = math.sqrt(accumulator["candidate_squared"])
        relative_l2 = difference_l2 / max(reference_l2, 1.0e-300)

        tolerance = tolerances[name]
        allowed_maximum = (
            tolerance.absolute
            + tolerance.relative
            * accumulator["maximum_reference_absolute"]
        )

        field_passed = (
            relative_l2 <= tolerance.relative
            and accumulator["maximum_absolute_difference"] <= allowed_maximum
        )

        all_passed = all_passed and field_passed

        report[name] = {
            "reference_l2": reference_l2,
            "candidate_l2": candidate_l2,
            "difference_l2": difference_l2,
            "relative_l2_difference": relative_l2,
            "maximum_absolute_difference": accumulator[
                "maximum_absolute_difference"
            ],
            "maximum_reference_absolute": accumulator[
                "maximum_reference_absolute"
            ],
            "allowed_maximum_difference": allowed_maximum,
            "passed": field_passed,
        }

    return all_passed, report


def read_residual_csv(output_dir, case_name):
    path = output_dir / "{}_distributed_residuals.csv".format(case_name)
    if not path.is_file():
        raise FileNotFoundError("Missing residual CSV: {}".format(path))

    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def compare_residual_histories(reference_rows, candidate_rows):
    common_count = min(len(reference_rows), len(candidate_rows))
    maximum_relative_difference = 0.0

    for index in range(common_count):
        reference_row = reference_rows[index]
        candidate_row = candidate_rows[index]

        for field in CSV_FIELDS:
            reference_value = float(reference_row[field])
            candidate_value = float(candidate_row[field])

            if not math.isfinite(reference_value) or not math.isfinite(
                candidate_value
            ):
                raise ValueError(
                    "Non-finite residual value at row {}, {}".format(
                        index + 1,
                        field,
                    )
                )

            scale = max(abs(reference_value), abs(candidate_value), 1.0e-300)
            maximum_relative_difference = max(
                maximum_relative_difference,
                abs(candidate_value - reference_value) / scale,
            )

    return {
        "reference_rows": len(reference_rows),
        "candidate_rows": len(candidate_rows),
        "common_rows": common_count,
        "maximum_scalar_relative_difference": maximum_relative_difference,
        "reference_cg_iterations": [
            int(row["cg_iterations"]) for row in reference_rows
        ],
        "candidate_cg_iterations": [
            int(row["cg_iterations"]) for row in candidate_rows
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_output", type=Path)
    parser.add_argument("candidate_output", type=Path)
    parser.add_argument("--case", default="blanket")
    parser.add_argument("--json-report", type=Path)
    args = parser.parse_args()

    reference_dir = args.reference_output.resolve()
    candidate_dir = args.candidate_output.resolve()

    if not reference_dir.is_dir():
        raise FileNotFoundError(
            "Reference output directory does not exist: {}".format(
                reference_dir
            )
        )

    if not candidate_dir.is_dir():
        raise FileNotFoundError(
            "Candidate output directory does not exist: {}".format(
                candidate_dir
            )
        )

    print("Loading reference final state")
    reference_iteration, reference_pieces, reference = load_owned_solution(
        reference_dir,
        args.case,
    )

    print("Loading candidate final state")
    candidate_iteration, candidate_pieces, candidate = load_owned_solution(
        candidate_dir,
        args.case,
    )

    if reference_iteration != candidate_iteration:
        raise ValueError(
            "Final output iterations differ: reference={}, candidate={}".format(
                reference_iteration,
                candidate_iteration,
            )
        )

    fields_passed, field_report = compare_fields(
        reference,
        candidate,
        DEFAULT_TOLERANCES,
    )

    residual_report = compare_residual_histories(
        read_residual_csv(reference_dir, args.case),
        read_residual_csv(candidate_dir, args.case),
    )

    report = {
        "case": args.case,
        "iteration": reference_iteration,
        "reference_piece_count": reference_pieces,
        "candidate_piece_count": candidate_pieces,
        "global_owned_node_count": len(reference),
        "fields_passed": fields_passed,
        "fields": field_report,
        "residual_history": residual_report,
    }

    print()
    print(
        "Compared {} unique owned global nodes at iteration {}".format(
            len(reference),
            reference_iteration,
        )
    )
    print(
        "Reference pieces: {}; candidate pieces: {}".format(
            reference_pieces,
            candidate_pieces,
        )
    )
    print()
    print(
        "{:>12} {:>16} {:>16} {:>16} {:>8}".format(
            "field",
            "relative L2",
            "maximum abs",
            "allowed max",
            "result",
        )
    )

    for name in FIELD_NAMES:
        item = field_report[name]
        print(
            "{:>12} {:16.8e} {:16.8e} {:16.8e} {:>8}".format(
                name,
                float(item["relative_l2_difference"]),
                float(item["maximum_absolute_difference"]),
                float(item["allowed_maximum_difference"]),
                "PASS" if item["passed"] else "FAIL",
            )
        )

    print()
    print(
        "Maximum residual-history scalar relative difference: {:.8e}".format(
            residual_report["maximum_scalar_relative_difference"]
        )
    )
    print(
        "CG iterations: {} versus {}".format(
            residual_report["reference_cg_iterations"],
            residual_report["candidate_cg_iterations"],
        )
    )

    if args.json_report:
        args.json_report.parent.mkdir(parents=True, exist_ok=True)
        args.json_report.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print("JSON report: {}".format(args.json_report))

    print()
    print(
        "DD-4D rank-count field comparison: {}".format(
            "PASS" if fields_passed else "FAIL"
        )
    )

    return 0 if fields_passed else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        raise SystemExit(2)
