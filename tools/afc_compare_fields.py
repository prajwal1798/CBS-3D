#!/usr/bin/env python3
"""Compare explicit and steady AFC PVTU fields on identical partitioned meshes.

Usage:
    python3 tools/afc_compare_fields.py EXPLICIT_DIR IMPLICIT_DIR --tol-k 1.0

Nonzero exit code means the fields fail temperature or velocity parity.
Both result arguments may be .pvtu files or directories containing the
latest .pvtu. The metrics include duplicated MPI ghost nodes; max norms
are unaffected by duplicates.
"""
import argparse
import sys

import numpy as np

from energy_balance import read_result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("explicit", help="Explicit reference PVTU or output directory")
    p.add_argument("implicit", help="Implicit candidate PVTU or output directory")
    p.add_argument("--tol-k", type=float, default=1.0,
                   help="Maximum allowed nodal |T_exp - T_imp| in K")
    p.add_argument("--velocity-tol", type=float, default=1e-12,
                   help="Maximum allowed nodal velocity error in m/s")
    a = p.parse_args()

    if a.tol_k < 0 or a.velocity_tol < 0:
        p.error("tolerances must be nonnegative")

    exp, path_e = read_result(a.explicit)
    imp, path_i = read_result(a.implicit)
    print("EXPLICIT:", path_e)
    print("IMPLICIT:", path_i)

    if len(exp) != len(imp):
        raise ValueError("Partition-piece counts differ; compare identical meshes")

    maximum_dt = -1.0
    maximum_du = 0.0
    sum_dt2 = 0.0
    count = 0
    hot_exp = (-np.inf, None)
    hot_imp = (-np.inf, None)
    maximum_point = None

    for piece_no, (e, i) in enumerate(zip(exp, imp)):
        xe, ce, ve, te, _ = e
        xi, ci, vi, ti, _ = i
        if xe.shape != xi.shape or ce.shape != ci.shape:
            raise ValueError(f"Partition {piece_no}: mesh shape differs")
        if not np.allclose(xe, xi, atol=1e-10, rtol=0):
            raise ValueError(f"Partition {piece_no}: node coordinates differ")
        if not np.array_equal(ce, ci):
            raise ValueError(f"Partition {piece_no}: cell connectivity differs")

        d = np.abs(te - ti)
        u = np.linalg.norm(ve - vi, axis=1)
        maximum_du = max(maximum_du, float(np.max(u)))
        sum_dt2 += float(np.dot(d, d))
        count += len(d)

        j = int(np.argmax(d))
        if float(d[j]) > maximum_dt:
            maximum_dt = float(d[j])
            maximum_point = (xe[j], float(te[j]), float(ti[j]))

        j = int(np.argmax(te))
        if float(te[j]) > hot_exp[0]:
            hot_exp = (float(te[j]), xe[j])

        j = int(np.argmax(ti))
        if float(ti[j]) > hot_imp[0]:
            hot_imp = (float(ti[j]), xi[j])

    print(f"MAX VELOCITY DIFFERENCE : {maximum_du:.9g} m/s")
    print(f"MAX TEMPERATURE ERROR   : {maximum_dt:.9g} K")
    print(f"RMS TEMPERATURE ERROR   : {(sum_dt2/count)**0.5:.9g} K")
    print(f"MAX ERROR AT            : {maximum_point[0]} m")
    print(f"T_EXPLICIT, T_IMPLICIT  : {maximum_point[1]:.9g}, {maximum_point[2]:.9g} K")
    print(f"EXPLICIT HOTSPOT        : {hot_exp[0]:.9g} K at {hot_exp[1]} m")
    print(f"IMPLICIT HOTSPOT        : {hot_imp[0]:.9g} K at {hot_imp[1]} m")
    passed = maximum_dt <= a.tol_k and maximum_du <= a.velocity_tol
    print(f"FIELD PARITY            : {'PASS' if passed else 'FAIL'} "
          f"(tolerance {a.tol_k:g} K, {a.velocity_tol:g} m/s)")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
