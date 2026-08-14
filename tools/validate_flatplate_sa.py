#!/usr/bin/env python3
"""
Flat-plate Spalart-Allmaras validation for CBS3D++_SI.

This script post-processes the VTU output of a turbulent flat-plate run and
compares it against the accepted reference behaviour for the case:

    1. Skin friction        cf(Re_x) against the Karman-Schoenherr and
                            1/7-power correlations.
    2. Law of the wall      u+(y+) against the viscous sublayer u+ = y+ and the
                            log law u+ = ln(y+)/kappa + B.
    3. Eddy viscosity       nu_t/nu profile shape and peak location.
    4. Wall resolution      y+ of the first cell, which must be below 1 for an
                            integrate-to-the-wall SA run.

The reference case is the NASA Turbulence Modeling Resource 2D zero-pressure
-gradient flat plate, run here as a thin three-dimensional slab.

Usage:
    python3 tools/validate_flatplate_sa.py \\
        --vtu out/flatplate_00500.vtu \\
        --rho 1.225 --mu 1.7894e-5 --uinf 34.6 \\
        --stations 0.2 0.5 0.97 \\
        --plate-origin 0.0 \\
        --out validation/

Requires: numpy, matplotlib, and either meshio or vtk.
"""

import argparse
import os
import sys

import numpy as np


KAPPA = 0.41
B_LOG = 5.0


def read_vtu(path):
    """Reads points and point-data arrays from a VTU file."""
    try:
        import meshio
    except ImportError:
        sys.exit(
            "meshio is required: pip install meshio\n"
            "(or adapt this reader to the vtk Python module)")

    mesh = meshio.read(path)
    return mesh.points, mesh.point_data


def find_array(point_data, candidates):
    """Returns the first present array among several possible field names."""
    for name in candidates:
        if name in point_data:
            return np.asarray(point_data[name])

    for name in candidates:
        for key in point_data:
            if key.lower() == name.lower():
                return np.asarray(point_data[key])

    return None


def extract_profile(points, values, x_station, tol):
    """Extracts a wall-normal profile at a streamwise station.

    Nodes within tol of x_station are collected and sorted by wall distance y.
    The spanwise direction is averaged over, which removes any residual
    three-dimensionality in the slab.
    """
    mask = np.abs(points[:, 0] - x_station) < tol

    if not np.any(mask):
        return None, None

    y = points[mask, 1]
    v = values[mask]

    order = np.argsort(y)
    y = y[order]
    v = v[order]

    # Average duplicate y values coming from different spanwise nodes.
    unique_y, inverse = np.unique(np.round(y, 12), return_inverse=True)
    averaged = np.zeros_like(unique_y)

    for i in range(len(unique_y)):
        averaged[i] = v[inverse == i].mean()

    return unique_y, averaged


def wall_shear_from_profile(y, u, mu):
    """Estimates the wall shear stress from the near-wall velocity profile.

    A one-sided three-point derivative on the graded mesh is used rather than a
    simple difference quotient, because the first cells are strongly stretched
    and a two-point estimate carries a first-order spacing error that shows up
    directly in cf.
    """
    if len(y) < 3:
        return None

    y0, y1, y2 = y[0], y[1], y[2]
    u0, u1, u2 = u[0], u[1], u[2]

    h1 = y1 - y0
    h2 = y2 - y0

    if h1 <= 0.0 or h2 <= 0.0 or abs(h2 - h1) < 1.0e-30:
        return None

    dudy = (
        -(h1 + h2) / (h1 * h2) * u0
        + h2 / (h1 * (h2 - h1)) * u1
        - h1 / (h2 * (h2 - h1)) * u2
    )

    return mu * dudy


def cf_karman_schoenherr(re_x):
    """Karman-Schoenherr skin friction, valid for turbulent boundary layers."""
    re_x = np.maximum(re_x, 1.0e3)
    # Solved iteratively: 1/sqrt(cf) = 4.13 log10(Re_x cf)
    cf = 0.0576 * re_x ** (-0.2)

    for _ in range(100):
        cf = 1.0 / (4.13 * np.log10(np.maximum(re_x * cf, 1.0e-12))) ** 2

    return cf


def cf_power_law(re_x):
    """1/7-power-law skin friction, a coarser but common reference."""
    return 0.0576 * np.maximum(re_x, 1.0e3) ** (-0.2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vtu", required=True, help="VTU file to analyse")
    parser.add_argument("--rho", type=float, required=True)
    parser.add_argument("--mu", type=float, required=True)
    parser.add_argument("--uinf", type=float, required=True)
    parser.add_argument("--plate-origin", type=float, default=0.0,
                        help="x of the leading edge, for Re_x")
    parser.add_argument("--stations", type=float, nargs="+",
                        default=[0.2, 0.5, 0.97])
    parser.add_argument("--tol", type=float, default=None,
                        help="streamwise tolerance for profile extraction")
    parser.add_argument("--out", default="validation")
    args = parser.parse_args()

    points, point_data = read_vtu(args.vtu)

    velocity = find_array(point_data, ["velocity", "U", "unkno", "Velocity"])
    nu_t = find_array(point_data, ["nu_t", "nut", "eddy_viscosity"])
    wall_distance = find_array(point_data, ["wall_distance", "d", "walldist"])

    if velocity is None:
        sys.exit("no velocity array found in the VTU point data")

    if velocity.ndim == 2:
        u = velocity[:, 0]
    else:
        u = velocity

    nu = args.mu / args.rho
    os.makedirs(args.out, exist_ok=True)

    if args.tol is None:
        span = points[:, 0].max() - points[:, 0].min()
        tol = 0.002 * span
    else:
        tol = args.tol

    print("Flat-plate Spalart-Allmaras validation")
    print("  file           %s" % args.vtu)
    print("  nodes          %d" % len(points))
    print("  nu             %.6e m^2/s" % nu)
    print("  U_inf          %.4f m/s" % args.uinf)
    print("  station tol    %.6e m" % tol)
    print("")

    rows = []
    profiles = {}

    for x_station in args.stations:
        y, u_profile = extract_profile(points, u, x_station, tol)

        if y is None:
            print("  station x=%.4f : no nodes found" % x_station)
            continue

        tau_w = wall_shear_from_profile(y, u_profile, args.mu)

        if tau_w is None or tau_w <= 0.0:
            print("  station x=%.4f : could not evaluate wall shear"
                  % x_station)
            continue

        u_tau = np.sqrt(tau_w / args.rho)
        cf = tau_w / (0.5 * args.rho * args.uinf ** 2)

        re_x = args.rho * args.uinf * (x_station - args.plate_origin) / args.mu

        cf_ks = cf_karman_schoenherr(re_x)
        cf_pl = cf_power_law(re_x)

        y_plus = y * u_tau / nu
        u_plus = u_profile / u_tau

        first_y_plus = y_plus[1] if len(y_plus) > 1 else y_plus[0]

        rows.append((x_station, re_x, cf, cf_ks, cf_pl, first_y_plus))
        profiles[x_station] = (y_plus, u_plus, y, u_profile)

        error_ks = 100.0 * (cf - cf_ks) / cf_ks

        print("  station x=%.4f" % x_station)
        print("    Re_x                 %.4e" % re_x)
        print("    cf (solver)          %.6e" % cf)
        print("    cf (Karman-Schoenherr) %.6e  (%+.2f %%)"
              % (cf_ks, error_ks))
        print("    cf (1/7 power law)   %.6e" % cf_pl)
        print("    u_tau                %.6f m/s" % u_tau)
        print("    first-cell y+        %.4f" % first_y_plus)

        if first_y_plus > 1.0:
            print("    WARNING: first-cell y+ exceeds 1; the near-wall SA")
            print("             solution is under-resolved and cf will be")
            print("             mesh dependent")

        # Log-law agreement in the region 30 < y+ < 300.
        band = (y_plus > 30.0) & (y_plus < 300.0)

        if np.count_nonzero(band) >= 3:
            log_law = np.log(y_plus[band]) / KAPPA + B_LOG
            deviation = np.max(np.abs(u_plus[band] - log_law))
            print("    max |u+ - log law|   %.3f  (30 < y+ < 300)" % deviation)
        else:
            print("    log-law band has too few points to assess")

        print("")

    if not rows:
        sys.exit("no usable stations; check --stations and --tol")

    # Text report.
    report = os.path.join(args.out, "flatplate_sa_report.csv")

    with open(report, "w") as handle:
        handle.write("x,Re_x,cf_solver,cf_karman_schoenherr,"
                     "cf_power_law,first_cell_yplus\n")

        for row in rows:
            handle.write("%.6f,%.6e,%.6e,%.6e,%.6e,%.4f\n" % row)

    print("  wrote %s" % report)

    # Plots.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available; skipping plots")
        return

    fig, axis = plt.subplots(figsize=(6.5, 4.5))

    for x_station, (y_plus, u_plus, _, _) in profiles.items():
        axis.semilogx(y_plus, u_plus, "o-", markersize=3,
                      label="CBS3D x=%.2f m" % x_station)

    reference = np.logspace(0, 3.5, 200)
    axis.semilogx(reference, reference, "k--", linewidth=1,
                  label=r"$u^+=y^+$")
    axis.semilogx(reference, np.log(reference) / KAPPA + B_LOG, "k:",
                  linewidth=1, label=r"$u^+=\ln(y^+)/\kappa+B$")

    axis.set_xlim(1.0, 3000.0)
    axis.set_ylim(0.0, 30.0)
    axis.set_xlabel(r"$y^+$")
    axis.set_ylabel(r"$u^+$")
    axis.set_title("Flat plate, law of the wall")
    axis.legend(fontsize=8)
    axis.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "law_of_the_wall.png"), dpi=150)

    fig, axis = plt.subplots(figsize=(6.5, 4.5))

    re_solver = np.array([r[1] for r in rows])
    cf_solver = np.array([r[2] for r in rows])

    re_reference = np.logspace(5, 7.5, 200)

    axis.loglog(re_solver, cf_solver, "o", markersize=6, label="CBS3D SA")
    axis.loglog(re_reference, cf_karman_schoenherr(re_reference), "k-",
                linewidth=1, label="Karman-Schoenherr")
    axis.loglog(re_reference, cf_power_law(re_reference), "k--",
                linewidth=1, label="1/7 power law")

    axis.set_xlabel(r"$Re_x$")
    axis.set_ylabel(r"$c_f$")
    axis.set_title("Flat plate, skin friction")
    axis.legend(fontsize=8)
    axis.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "skin_friction.png"), dpi=150)

    print("  wrote plots to %s" % args.out)

    if nu_t is not None:
        peak = np.max(nu_t) / nu
        print("")
        print("  peak nu_t/nu   %.1f" % peak)
        print("  (a converged flat-plate SA solution reaches order 100-1000")
        print("   depending on Re; a value near zero means the SA source terms")
        print("   are not firing, and a runaway value means the destruction")
        print("   term or the wall distance is wrong)")

    if wall_distance is not None:
        print("  wall distance  min %.4e  max %.4e"
              % (np.min(wall_distance), np.max(wall_distance)))


if __name__ == "__main__":
    main()
