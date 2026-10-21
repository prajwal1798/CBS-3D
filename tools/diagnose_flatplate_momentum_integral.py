#!/usr/bin/env python3
"""
diagnose_flatplate_momentum_integral.py

Distributed CBS3D flat-plate boundary-layer integral audit.

Purpose
-------
Uses the existing, production CBS3D distributed TMR reader so that the
diagnostic sees exactly the same rank-local mesh, ownership metadata, wall
faces, and TET4 velocity gradients as analyze_distributed_tmr.py.

For each exact streamwise mesh plane in the requested x-range it reports:
  - U_e from the outermost wall-normal profile levels
  - delta*, theta, delta99, H12
  - Re_theta and 0.036 Re_x^0.8
  - wall Cf from:
      (1) TET4 normal velocity gradient
      (2) full deviatoric TET4 stress traction
      (3) independent 3-point profile wall derivative
  - ZPG momentum-integral residual:
        Cf/2 - dtheta/dx
    using fixed U_inf
  - generalized incompressible momentum-integral residual:
        Cf/2 - [dtheta/dx + (delta*+2 theta)/U_e dU_e/dx]
    using locally measured U_e.

No third-party Python package is required.  Compatible with Python 3.6+.
"""

from __future__ import print_function

import argparse
import csv
import importlib.util
import json
import math
import re
import statistics
from pathlib import Path


def load_base(repo_root):
    path = (repo_root / "examples" / "Flat_Plate_Turbulent" /
            "analyze_distributed_tmr.py")
    if not path.is_file():
        raise SystemExit("missing base analyzer: {}".format(path))
    spec = importlib.util.spec_from_file_location(
        "cbs_flatplate_analyzer", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("partition_root", type=Path)
    p.add_argument("output_dir", type=Path)
    p.add_argument("--repo-root", type=Path,
                   default=Path(__file__).resolve().parents[1])
    p.add_argument("--u-inf", type=float, default=1.0)
    p.add_argument("--re-x1", type=float, default=5.0e6)
    p.add_argument("--x-min", type=float, default=0.05)
    p.add_argument("--x-max", type=float, default=1.90)
    p.add_argument("--outer-levels", type=int, default=8,
                   help="highest wall-normal levels averaged for local U_e")
    p.add_argument("--station", type=float, default=0.97008)
    p.add_argument("--csv", type=Path, default=None)
    p.add_argument("--json", type=Path, default=None)
    return p.parse_args()


def round_x(x):
    return round(float(x), 12)


def round_y(y):
    return round(float(y), 14)


def trapz_profile(y, f):
    total = 0.0
    for i in range(1, len(y)):
        total += 0.5 * (f[i - 1] + f[i]) * (y[i] - y[i - 1])
    return total


def station_integrals(y, u, ue):
    if len(y) < 3:
        raise ValueError("profile has fewer than 3 points")
    if not (ue > 0.0) or not math.isfinite(ue):
        raise ValueError("invalid edge velocity")

    r = [v / ue for v in u]
    dstar = trapz_profile(y, [1.0 - q for q in r])
    theta = trapz_profile(y, [q * (1.0 - q) for q in r])

    delta99 = float("nan")
    target = 0.99
    for i in range(1, len(y)):
        if r[i - 1] < target <= r[i]:
            den = r[i] - r[i - 1]
            if abs(den) > 1.0e-300:
                f = (target - r[i - 1]) / den
                delta99 = y[i - 1] + f * (y[i] - y[i - 1])
            else:
                delta99 = y[i]
            break

    return {
        "delta_star": dstar,
        "theta": theta,
        "delta99": delta99,
        "H12": dstar / theta if theta > 0.0 else float("nan"),
    }


def wall_dudy_three_point(y, u):
    if len(y) < 3:
        raise ValueError("need at least 3 profile points")
    if abs(y[0]) > 1.0e-12:
        raise ValueError("profile has no y=0 wall node")

    h1 = y[1] - y[0]
    h2 = y[2] - y[0]
    if not (h1 > 0.0 and h2 > h1):
        raise ValueError("invalid first three wall-normal spacings")

    return (
        -(h1 + h2) / (h1 * h2) * u[0]
        + h2 / (h1 * (h2 - h1)) * u[1]
        - h1 / (h2 * (h2 - h1)) * u[2]
    )


def derivative_three_point(xs, fs):
    """Second-order 3-point derivative on an arbitrary monotone x grid."""
    n = len(xs)
    if n < 3:
        raise ValueError("need >=3 stations for derivative")

    out = [0.0] * n

    def lagrange_derivative_at(indices, xeval):
        coeff = []
        for jj in indices:
            xj = xs[jj]
            c = 0.0
            for mm in indices:
                if mm == jj:
                    continue
                term = 1.0 / (xj - xs[mm])
                for kk in indices:
                    if kk == jj or kk == mm:
                        continue
                    term *= (xeval - xs[kk]) / (xj - xs[kk])
                c += term
            coeff.append(c)
        return sum(c * fs[j] for c, j in zip(coeff, indices))

    out[0] = lagrange_derivative_at([0, 1, 2], xs[0])
    for i in range(1, n - 1):
        out[i] = lagrange_derivative_at([i - 1, i, i + 1], xs[i])
    out[-1] = lagrange_derivative_at([n - 3, n - 2, n - 1], xs[-1])
    return out


def linear_interp(xvals, yvals, x):
    if not xvals:
        return None
    if x < xvals[0] or x > xvals[-1]:
        return None
    if x == xvals[0]:
        return yvals[0]
    if x == xvals[-1]:
        return yvals[-1]
    lo = 0
    hi = len(xvals) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xvals[mid] <= x:
            lo = mid
        else:
            hi = mid
    x0, x1 = xvals[lo], xvals[hi]
    if x1 == x0:
        return yvals[lo]
    f = (x - x0) / (x1 - x0)
    return yvals[lo] + f * (yvals[hi] - yvals[lo])


def vector_norm(v):
    return math.sqrt(sum(q * q for q in v))


def vector_dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def vector_tangent(v, n):
    vn = vector_dot(v, n)
    return [v[i] - vn * n[i] for i in range(3)]


def group_wall_samples(samples):
    grouped = {}
    for x, area, simple_shear, dev_shear, simple_eff_shear, dev_eff_shear in samples:
        key = round_x(x)
        g = grouped.setdefault(key, [0.0, 0.0, 0.0, 0.0, 0.0])
        g[0] += area
        g[1] += area * simple_shear
        g[2] += area * dev_shear
        g[3] += area * simple_eff_shear
        g[4] += area * dev_eff_shear

    xs = sorted(grouped)
    simple = []
    dev = []
    simple_eff = []
    dev_eff = []
    for x in xs:
        area, ssum, dsum, sesum, desum = grouped[x]
        simple.append(ssum / area)
        dev.append(dsum / area)
        simple_eff.append(sesum / area)
        dev_eff.append(desum / area)
    return xs, simple, dev, simple_eff, dev_eff


def rel_error(lhs, rhs):
    return abs(lhs - rhs) / max(abs(lhs), 1.0e-300)


def median_finite(values):
    finite = [v for v in values if math.isfinite(v)]
    return statistics.median(finite) if finite else float("nan")


def mean(values):
    return sum(values) / float(len(values)) if values else float("nan")


def main():
    args = parse_args()
    base = load_base(args.repo_root)

    if not (args.u_inf > 0.0):
        raise SystemExit("--u-inf must be positive")
    if not (args.re_x1 > 0.0):
        raise SystemExit("--re-x1 must be positive")
    if args.outer_levels < 2:
        raise SystemExit("--outer-levels must be >=2")

    nu = args.u_inf / args.re_x1

    ranks = sorted(
        p for p in args.partition_root.glob("rank_[0-9][0-9][0-9][0-9]")
        if p.is_dir()
    )
    if len(ranks) != 40:
        raise SystemExit(
            "expected 40 partition ranks, found {}".format(len(ranks)))

    rx = re.compile(r"flatplate_step_(\d{8})_rank_(\d{4})\.vtu$")
    by_step = {}
    for f in args.output_dir.glob("flatplate_step_*_rank_*.vtu"):
        m = rx.match(f.name)
        if m:
            by_step.setdefault(int(m.group(1)), {})[int(m.group(2))] = f

    good = [step for step, pieces in by_step.items() if len(pieces) == 40]
    if not good:
        raise SystemExit("no complete 40-piece distributed snapshot")
    step = max(good)

    owned = {}
    wall_samples = []

    for irank, rank in enumerate(ranks):
        tets, mesh_points, bfaces, mapping = base.read_partition(rank)
        points, velocity, pd, cd = base.read_piece(by_step[step][irank])

        if len(points) != len(mesh_points):
            raise RuntimeError(
                "mesh/VTU point mismatch on rank {}".format(irank))

        gids = [int(round(v)) for v in pd["global_node_id"]]
        is_owned = [int(round(v)) for v in pd["is_owned"]]

        for i, gid in enumerate(gids):
            if is_owned[i] == 1:
                owned[gid] = (points[i], velocity[i])

        for face, parent, rawbc in bfaces:
            if mapping.get(rawbc) != 530:
                continue

            pts = [mesh_points[n - 1] for n in face]
            xc = sum(p[0] for p in pts) / 3.0
            if xc <= 0.0:
                continue

            tet = tets[parent - 1]
            face_set = set(face)
            opp = [n for n in tet if n not in face_set]
            if len(opp) != 1:
                raise RuntimeError(
                    "wall face/parent mismatch rank {} parent {}".format(
                        irank, parent))

            p0, p1, p2 = pts
            normal = base.cross(base.sub(p1, p0), base.sub(p2, p0))
            centroid = (
                sum(p[0] for p in pts) / 3.0,
                sum(p[1] for p in pts) / 3.0,
                sum(p[2] for p in pts) / 3.0,
            )
            opposite = mesh_points[opp[0] - 1]
            if base.dot(normal, base.sub(opposite, centroid)) > 0.0:
                normal = (-normal[0], -normal[1], -normal[2])

            nn = base.norm(normal)
            if not (nn > 0.0):
                raise RuntimeError("zero wall-face normal")
            nunit = [normal[i] / nn for i in range(3)]
            area = 0.5 * nn

            gu = base.velocity_gradient(mesh_points, tet, velocity)

            # Existing TMR analyzer's one-sided normal-gradient traction.
            simple_traction = [
                nu * vector_dot(gu[i], nunit) for i in range(3)
            ]
            simple_shear = vector_norm(
                vector_tangent(simple_traction, nunit))

            # Full deviatoric kinematic stress used by the RANS momentum model.
            divu = gu[0][0] + gu[1][1] + gu[2][2]
            tau = [[0.0] * 3 for _ in range(3)]
            for i in range(3):
                for j in range(3):
                    tau[i][j] = nu * (
                        gu[i][j] + gu[j][i]
                        - ((2.0 / 3.0) * divu if i == j else 0.0)
                    )
            dev_traction = [
                sum(tau[i][j] * nunit[j] for j in range(3))
                for i in range(3)
            ]
            dev_shear = vector_norm(vector_tangent(dev_traction, nunit))

            # The production momentum operator uses element effective viscosity.
            # Report that discrete wall-adjacent stress separately from the
            # physical molecular wall stress.  This catches eddy-viscosity
            # leakage into the first TET4 as a distinct numerical effect.
            mu_e = cd.get("mu_e", [])
            mu_t_e = cd.get("mu_t_e", [])
            eff_ratio = 1.0
            if mu_e and mu_t_e:
                me = mu_e[parent - 1]
                mte = mu_t_e[parent - 1]
                if me > 0.0 and math.isfinite(me) and math.isfinite(mte):
                    eff_ratio = (me + mte) / me

            simple_eff_shear = simple_shear * eff_ratio
            dev_eff_shear = dev_shear * eff_ratio

            wall_samples.append(
                (xc, area, simple_shear, dev_shear,
                 simple_eff_shear, dev_eff_shear))

    if len(owned) < 1:
        raise RuntimeError("no owned nodes recovered")
    if len(wall_samples) < 1:
        raise RuntimeError("no flat-plate wall faces recovered")

    # Exact streamwise mesh planes from owned nodes. Average spanwise copies
    # only after grouping by exact x and y.
    planes = {}
    for coord, vel in owned.values():
        x, y = coord[0], coord[1]
        if x < args.x_min - 1.0e-12 or x > args.x_max + 1.0e-12:
            continue
        xk = round_x(x)
        yk = round_y(y)
        g = planes.setdefault(xk, {}).setdefault(yk, [0.0, 0])
        g[0] += vel[0]
        g[1] += 1

    wall_x, wall_simple, wall_dev, wall_simple_eff, wall_dev_eff = (
        group_wall_samples(wall_samples)
    )

    rows = []
    for x in sorted(planes):
        ygroups = planes[x]
        prof = []
        for y in sorted(ygroups):
            usum, count = ygroups[y]
            prof.append((float(y), usum / float(count)))

        if len(prof) < max(5, args.outer_levels):
            continue
        if abs(prof[0][0]) > 1.0e-12:
            continue

        y = [p[0] for p in prof]
        u = [p[1] for p in prof]

        outer = u[-min(args.outer_levels, len(u)):]
        ue = mean(outer)
        if not (ue > 0.0 and math.isfinite(ue)):
            continue

        fixed = station_integrals(y, u, args.u_inf)
        local = station_integrals(y, u, ue)

        try:
            dudy = wall_dudy_three_point(y, u)
        except ValueError:
            continue

        shear_simple = linear_interp(wall_x, wall_simple, x)
        shear_dev = linear_interp(wall_x, wall_dev, x)
        shear_simple_eff = linear_interp(wall_x, wall_simple_eff, x)
        shear_dev_eff = linear_interp(wall_x, wall_dev_eff, x)
        if (shear_simple is None or shear_dev is None or
                shear_simple_eff is None or shear_dev_eff is None):
            continue

        cf_profile_inf = 2.0 * nu * abs(dudy) / (args.u_inf ** 2)
        cf_simple_inf = 2.0 * shear_simple / (args.u_inf ** 2)
        cf_dev_inf = 2.0 * shear_dev / (args.u_inf ** 2)
        cf_simple_eff_inf = 2.0 * shear_simple_eff / (args.u_inf ** 2)
        cf_dev_eff_inf = 2.0 * shear_dev_eff / (args.u_inf ** 2)

        cf_profile_edge = 2.0 * nu * abs(dudy) / (ue ** 2)
        cf_simple_edge = 2.0 * shear_simple / (ue ** 2)
        cf_dev_edge = 2.0 * shear_dev / (ue ** 2)
        cf_simple_eff_edge = 2.0 * shear_simple_eff / (ue ** 2)
        cf_dev_eff_edge = 2.0 * shear_dev_eff / (ue ** 2)

        rex = args.u_inf * x / nu
        re_theta = args.u_inf * fixed["theta"] / nu
        re_theta_ref = 0.036 * (rex ** 0.8)

        rows.append({
            "x": x,
            "ue": ue,

            "theta_inf": fixed["theta"],
            "delta_star_inf": fixed["delta_star"],
            "delta99_inf": fixed["delta99"],
            "H12_inf": fixed["H12"],

            "theta_edge": local["theta"],
            "delta_star_edge": local["delta_star"],
            "delta99_edge": local["delta99"],
            "H12_edge": local["H12"],

            "cf_profile_inf": cf_profile_inf,
            "cf_tet_grad_inf": cf_simple_inf,
            "cf_tet_dev_inf": cf_dev_inf,
            "cf_tet_grad_effective_inf": cf_simple_eff_inf,
            "cf_tet_dev_effective_inf": cf_dev_eff_inf,

            "cf_profile_edge": cf_profile_edge,
            "cf_tet_grad_edge": cf_simple_edge,
            "cf_tet_dev_edge": cf_dev_edge,
            "cf_tet_grad_effective_edge": cf_simple_eff_edge,
            "cf_tet_dev_effective_edge": cf_dev_eff_edge,

            "Re_x": rex,
            "Re_theta": re_theta,
            "Re_theta_ref": re_theta_ref,
            "Re_theta_ratio": (
                re_theta / re_theta_ref
                if re_theta_ref > 0.0 else float("nan")
            ),
        })

    if len(rows) < 5:
        raise RuntimeError(
            "too few usable exact streamwise planes ({})".format(len(rows)))

    xs = [r["x"] for r in rows]
    dtheta_inf = derivative_three_point(
        xs, [r["theta_inf"] for r in rows])
    dtheta_edge = derivative_three_point(
        xs, [r["theta_edge"] for r in rows])
    due_dx = derivative_three_point(
        xs, [r["ue"] for r in rows])

    for i, r in enumerate(rows):
        lhs_zpg = 0.5 * r["cf_tet_dev_inf"]
        rhs_zpg = dtheta_inf[i]

        lhs_gen = 0.5 * r["cf_tet_dev_edge"]
        pressure_term = (
            (r["delta_star_edge"] + 2.0 * r["theta_edge"])
            / r["ue"] * due_dx[i]
        )
        rhs_gen = dtheta_edge[i] + pressure_term

        r["dtheta_inf_dx"] = dtheta_inf[i]
        r["dtheta_edge_dx"] = dtheta_edge[i]
        r["dUe_dx"] = due_dx[i]
        r["outer_velocity_term"] = pressure_term

        r["zpg_lhs_cf_over_2"] = lhs_zpg
        r["zpg_rhs_dtheta_dx"] = rhs_zpg
        r["zpg_rel_error"] = rel_error(lhs_zpg, rhs_zpg)

        r["general_lhs_cf_over_2"] = lhs_gen
        r["general_rhs"] = rhs_gen
        r["general_rel_error"] = rel_error(lhs_gen, rhs_gen)

        r["cf_tet_grad_vs_dev_rel"] = rel_error(
            r["cf_tet_dev_inf"], r["cf_tet_grad_inf"])
        r["cf_profile_vs_dev_rel"] = rel_error(
            r["cf_tet_dev_inf"], r["cf_profile_inf"])
        r["cf_effective_vs_molecular_rel"] = rel_error(
            r["cf_tet_dev_inf"], r["cf_tet_dev_effective_inf"])

    interior = rows[1:-1]

    summary = {
        "iteration": step,
        "owned_global_nodes": len(owned),
        "streamwise_plane_count": len(rows),
        "x_range": [rows[0]["x"], rows[-1]["x"]],
        "nu": nu,
        "u_inf": args.u_inf,

        "median_zpg_momentum_integral_rel_error_interior":
            median_finite([r["zpg_rel_error"] for r in interior]),

        "median_generalized_momentum_integral_rel_error_interior":
            median_finite([r["general_rel_error"] for r in interior]),

        "median_cf_tet_grad_vs_deviatoric_rel_difference":
            median_finite([r["cf_tet_grad_vs_dev_rel"] for r in rows]),

        "median_cf_profile_vs_deviatoric_rel_difference":
            median_finite([r["cf_profile_vs_dev_rel"] for r in rows]),

        "median_wall_effective_vs_molecular_cf_rel_difference":
            median_finite([r["cf_effective_vs_molecular_rel"] for r in rows]),

        "mean_Ue": mean([r["ue"] for r in rows]),
        "min_Ue": min(r["ue"] for r in rows),
        "max_Ue": max(r["ue"] for r in rows),
    }

    # Selected stations nearest the requested engineering locations.
    targets = [0.10, 0.30, 0.50, args.station, 1.90]
    selected = []
    used = set()
    for target in targets:
        idx = min(range(len(rows)), key=lambda i: abs(rows[i]["x"] - target))
        if idx in used:
            continue
        used.add(idx)
        r = rows[idx]
        selected.append({
            "target_x": target,
            "x": r["x"],
            "Re_x": r["Re_x"],
            "Re_theta": r["Re_theta"],
            "Re_theta_ref": r["Re_theta_ref"],
            "Re_theta_ratio": r["Re_theta_ratio"],
            "theta_inf": r["theta_inf"],
            "delta_star_inf": r["delta_star_inf"],
            "delta99_inf": r["delta99_inf"],
            "H12_inf": r["H12_inf"],
            "Ue": r["ue"],
            "cf_tet_dev_inf": r["cf_tet_dev_inf"],
            "cf_profile_inf": r["cf_profile_inf"],
            "zpg_rel_error": r["zpg_rel_error"],
            "general_rel_error": r["general_rel_error"],
        })
    summary["selected_stations"] = selected

    print("=" * 78)
    print("CBS3D DISTRIBUTED FLAT-PLATE MOMENTUM-INTEGRAL AUDIT")
    print("=" * 78)
    print("iteration                         : {}".format(step))
    print("exact streamwise planes           : {}".format(len(rows)))
    print("x range                           : {:.9g} .. {:.9g}".format(
        rows[0]["x"], rows[-1]["x"]))
    print("Ue min / mean / max               : {:.9g} / {:.9g} / {:.9g}".format(
        summary["min_Ue"], summary["mean_Ue"], summary["max_Ue"]))
    print()
    print("median ZPG residual               : {:.3f} %".format(
        100.0 * summary[
            "median_zpg_momentum_integral_rel_error_interior"]))
    print("median generalized residual       : {:.3f} %".format(
        100.0 * summary[
            "median_generalized_momentum_integral_rel_error_interior"]))
    print()
    print("median Cf TET-grad vs deviatoric  : {:.3f} %".format(
        100.0 * summary[
            "median_cf_tet_grad_vs_deviatoric_rel_difference"]))
    print("median Cf profile vs deviatoric   : {:.3f} %".format(
        100.0 * summary[
            "median_cf_profile_vs_deviatoric_rel_difference"]))
    print("median wall effective vs molecular: {:.3f} %".format(
        100.0 * summary[
            "median_wall_effective_vs_molecular_cf_rel_difference"]))

    print()
    print("SELECTED STATIONS")
    print(" x        Re_theta   ref       ratio    Cf(dev)    H12      "
          "ZPGerr%   GENerr%   Ue")
    for r in selected:
        print(
            "{:7.4f}  {:9.1f} {:9.1f} {:7.3f}  {:9.6f}  {:7.4f}  "
            "{:8.2f}  {:8.2f}  {:9.6f}".format(
                r["x"], r["Re_theta"], r["Re_theta_ref"],
                r["Re_theta_ratio"], r["cf_tet_dev_inf"], r["H12_inf"],
                100.0 * r["zpg_rel_error"],
                100.0 * r["general_rel_error"], r["Ue"])
        )

    zpg = summary["median_zpg_momentum_integral_rel_error_interior"]
    gen = summary["median_generalized_momentum_integral_rel_error_interior"]
    print()
    if gen < 0.05:
        verdict = (
            "GENERALIZED MOMENTUM INTEGRAL IS INTERNALLY CONSISTENT (<5%). "
            "Wall shear and integral momentum-thickness growth agree."
        )
    elif gen < 0.20:
        verdict = (
            "GENERALIZED MOMENTUM INTEGRAL IS MARGINAL (5-20%). "
            "Inspect stationwise residuals and x-resolution before changing physics."
        )
    else:
        verdict = (
            "GENERALIZED MOMENTUM INTEGRAL IS INCONSISTENT (>20%). "
            "Wall shear, momentum development, convergence, or boundary forcing "
            "must be resolved before blaming SA."
        )
    print("VERDICT: {}".format(verdict))

    if abs(zpg - gen) > 0.05:
        print(
            "NOTE: ZPG and generalized residuals differ materially; "
            "local Ue variation is not negligible."
        )

    if args.csv is None:
        csv_path = args.output_dir.parent / "momentum_integral_audit.csv"
    else:
        csv_path = args.csv
    if args.json is None:
        json_path = args.output_dir.parent / "momentum_integral_summary.json"
    else:
        json_path = args.json

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "x", "ue", "Re_x", "Re_theta", "Re_theta_ref", "Re_theta_ratio",
        "theta_inf", "delta_star_inf", "delta99_inf", "H12_inf",
        "theta_edge", "delta_star_edge", "delta99_edge", "H12_edge",
        "cf_profile_inf", "cf_tet_grad_inf", "cf_tet_dev_inf",
        "cf_tet_grad_effective_inf", "cf_tet_dev_effective_inf",
        "cf_profile_edge", "cf_tet_grad_edge", "cf_tet_dev_edge",
        "cf_tet_grad_effective_edge", "cf_tet_dev_effective_edge",
        "dtheta_inf_dx", "dtheta_edge_dx", "dUe_dx",
        "outer_velocity_term",
        "zpg_lhs_cf_over_2", "zpg_rhs_dtheta_dx", "zpg_rel_error",
        "general_lhs_cf_over_2", "general_rhs", "general_rel_error",
        "cf_tet_grad_vs_dev_rel", "cf_profile_vs_dev_rel",
        "cf_effective_vs_molecular_rel",
    ]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in fields})

    with json_path.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print()
    print("CSV  = {}".format(csv_path))
    print("JSON = {}".format(json_path))


if __name__ == "__main__":
    main()
