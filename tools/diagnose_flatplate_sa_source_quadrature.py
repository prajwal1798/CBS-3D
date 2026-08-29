#!/usr/bin/env python3
"""Frozen-field audit of the nonlinear SA source integration on the CBS3D flat plate.

This script does not run the CFD solver.  It reads one complete distributed VTU
snapshot plus the immutable rank-local partitions and evaluates the SA
production/destruction source two ways on the *same* frozen field:

1. current Liu production path:
       closure coefficients frozen at element-average q and d,
       remaining q_h integrated with the consistent P1 mass/load;

2. positive four-point tetrahedral quadrature:
       q_h and d_h interpolated to four quadrature points and the nonlinear
       closure evaluated pointwise.

Everything else is identical: SA-noft2, standard |curl u| vorticity, the same
molecular viscosity, the same r clamp, and the same P1 velocity/q fields.

The result localises whether the element-frozen nonlinear source closure can
plausibly explain the excessive outer-layer nu_tilde / mu_t seen in the NASA
flat-plate benchmark.
"""

from __future__ import print_function

import argparse
import importlib.util
import json
import math
from pathlib import Path

CB1 = 0.1355
CB2 = 0.622
SIGMA = 2.0 / 3.0
KAPPA = 0.41
CW2 = 0.3
CW3 = 2.0
CV1 = 7.1
CW1 = CB1 / (KAPPA * KAPPA) + (1.0 + CB2) / SIGMA

TET_A = 0.58541019662496845446
TET_B = 0.13819660112501051518
TET_W = 0.25

BINS = (
    (0.0, 1.0e-4),
    (1.0e-4, 5.0e-4),
    (5.0e-4, 1.0e-3),
    (1.0e-3, 5.0e-3),
    (5.0e-3, 1.0e-2),
    (1.0e-2, 2.0e-2),
    (2.0e-2, 5.0e-2),
    (5.0e-2, float("inf")),
)


def load_base_analyzer(repo_root):
    path = repo_root / "examples" / "Flat_Plate_Turbulent" / "analyze_distributed_tmr.py"
    if not path.is_file():
        raise SystemExit("missing base analyzer: {}".format(path))
    spec = importlib.util.spec_from_file_location("cbs_flatplate_analyzer", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("partition_root", type=Path)
    p.add_argument("output_dir", type=Path)
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--station", type=float, default=0.97008)
    p.add_argument("--station-half-width", type=float, default=0.03)
    p.add_argument("--min-wall-distance", type=float, default=1.0e-14)
    p.add_argument("--min-stilde", type=float, default=1.0e-14)
    p.add_argument("--top", type=int, default=12)
    return p.parse_args()


def fv1(chi):
    if not (chi > 0.0) or not math.isfinite(chi):
        return 0.0
    if chi > 1.0e6:
        r = CV1 / chi
        return 1.0 / (1.0 + r * r * r)
    chi3 = chi * chi * chi
    return chi3 / (chi3 + CV1 * CV1 * CV1)


def fv2(chi):
    f1 = fv1(chi)
    den = 1.0 + chi * f1
    return 0.0 if abs(den) <= 1.0e-30 else 1.0 - chi / den


def fw_from_r(r):
    r = max(0.0, min(10.0, r))
    r6 = r ** 6
    g = r + CW2 * (r6 - r)
    g6 = g ** 6
    cw36 = CW3 ** 6
    return g * ((1.0 + cw36) / (g6 + cw36)) ** (1.0 / 6.0)


def closure(q, d, omega, nu, min_d, min_stilde):
    d = max(d, min_d)
    q = max(q, 0.0)
    chi = q / nu
    sbar = 0.0
    if q > 0.0:
        sbar = q * fv2(chi) / (KAPPA * KAPPA * d * d)
    stilde = max(omega + sbar, min_stilde)
    den = stilde * KAPPA * KAPPA * d * d
    rval = 10.0 if den <= 1.0e-30 else max(0.0, min(10.0, q / den))
    fw = fw_from_r(rval)
    p = CB1 * stilde * q
    cd = CW1 * fw / (d * d)  # SA-noft2
    dterm = cd * q * q
    return p, dterm, stilde, rval, fw


def vorticity(base, points, tet, velocity):
    gu = base.velocity_gradient(points, tet, velocity)
    wx = gu[2][1] - gu[1][2]
    wy = gu[0][2] - gu[2][0]
    wz = gu[1][0] - gu[0][1]
    return math.sqrt(wx * wx + wy * wy + wz * wz)


def tet_volume(base, points, tet):
    p0, p1, p2, p3 = (points[n - 1] for n in tet)
    a = base.sub(p1, p0)
    b = base.sub(p2, p0)
    c = base.sub(p3, p0)
    return abs(base.dot(a, base.cross(b, c))) / 6.0


def consistent_q_load(volume, q, a):
    return volume * (sum(q) + q[a]) / 20.0


def current_liu_source(volume, q, d, omega, nu, min_d, min_stilde):
    qbar = 0.25 * sum(q)
    dbar = max(0.25 * sum(d), min_d)
    pbar, dbar_term, stilde, rval, fw = closure(
        qbar, dbar, omega, nu, min_d, min_stilde
    )
    pcoef = 0.0 if qbar <= 0.0 else pbar / qbar
    dcoef = 0.0 if qbar <= 0.0 else dbar_term / (qbar * qbar)

    pv = [0.0] * 4
    dv = [0.0] * 4
    for a in range(4):
        qload = consistent_q_load(volume, q, a)
        pv[a] = pcoef * qload
        dv[a] = dcoef * qbar * qload

    return {
        "p_vec": pv,
        "d_vec": dv,
        "p_total": sum(pv),
        "d_total": sum(dv),
        "qbar": qbar,
        "dbar": dbar,
        "stilde": stilde,
        "r": rval,
        "fw": fw,
    }


def quadrature_source(volume, q, d, omega, nu, min_d, min_stilde):
    pv = [0.0] * 4
    dv = [0.0] * 4
    qsets = (
        (TET_A, TET_B, TET_B, TET_B),
        (TET_B, TET_A, TET_B, TET_B),
        (TET_B, TET_B, TET_A, TET_B),
        (TET_B, TET_B, TET_B, TET_A),
    )

    for n in qsets:
        qq = sum(n[a] * q[a] for a in range(4))
        dd = sum(n[a] * d[a] for a in range(4))
        pp, dst, _, _, _ = closure(qq, dd, omega, nu, min_d, min_stilde)
        weight = volume * TET_W
        for a in range(4):
            pv[a] += weight * n[a] * pp
            dv[a] += weight * n[a] * dst

    return {
        "p_vec": pv,
        "d_vec": dv,
        "p_total": sum(pv),
        "d_total": sum(dv),
    }


def new_stats(lo, hi):
    return {
        "d_lo": lo,
        "d_hi": hi,
        "count": 0,
        "volume": 0.0,
        "p_current": 0.0,
        "d_current": 0.0,
        "p_quad": 0.0,
        "d_quad": 0.0,
        "err_sum": 0.0,
        "err_max": 0.0,
    }


def safe_ratio(a, b):
    if abs(b) <= 1.0e-300:
        return float("nan")
    return a / b


def main():
    args = parse_args()
    base = load_base_analyzer(args.repo_root)

    ranks = sorted(
        p for p in args.partition_root.glob("rank_[0-9][0-9][0-9][0-9]")
        if p.is_dir()
    )
    if len(ranks) != 40:
        raise SystemExit("expected 40 partition ranks, found {}".format(len(ranks)))

    rx = base.re.compile(r"flatplate_step_(\d{8})_rank_(\d{4})\.vtu$")
    by_step = {}
    for f in args.output_dir.glob("flatplate_step_*_rank_*.vtu"):
        m = rx.match(f.name)
        if m:
            by_step.setdefault(int(m.group(1)), {})[int(m.group(2))] = f
    good = [s for s, pieces in by_step.items() if len(pieces) == 40]
    if not good:
        raise SystemExit("no complete 40-piece distributed snapshot")
    step = max(good)

    bins = [new_stats(lo, hi) for lo, hi in BINS]
    total = new_stats(0.0, float("inf"))
    worst = []
    owned_profile = {}
    station_element_peak = 0.0
    station_element_peak_y = float("nan")

    for irank, rank in enumerate(ranks):
        tets, mesh_points, _, _ = base.read_partition(rank)
        points, velocity, pd, cd = base.read_piece(by_step[step][irank])
        if len(points) != len(mesh_points):
            raise RuntimeError("mesh/VTU point mismatch rank {}".format(irank))

        qfield = pd["nu_tilde"]
        dfield = pd["wall_distance"]
        ratios = pd.get("mu_t_over_mu", [])
        gids = [int(round(x)) for x in pd["global_node_id"]]
        owned = [int(round(x)) for x in pd["is_owned"]]

        for i, gid in enumerate(gids):
            if owned[i] == 1:
                owned_profile[gid] = (
                    points[i],
                    qfield[i],
                    ratios[i] if ratios else float("nan"),
                )

        mu_e = cd["mu_e"]
        rho_e = cd["rho_e"]
        mu_t_e = cd["mu_t_e"]

        for ie, tet in enumerate(tets):
            xc = sum(points[n - 1][0] for n in tet) * 0.25
            if abs(xc - args.station) > args.station_half_width:
                continue

            yc = sum(points[n - 1][1] for n in tet) * 0.25
            q = [qfield[n - 1] for n in tet]
            d = [dfield[n - 1] for n in tet]
            dc = 0.25 * sum(d)
            volume = tet_volume(base, points, tet)
            if not (mu_e[ie] > 0.0 and rho_e[ie] > 0.0):
                continue
            nu = mu_e[ie] / rho_e[ie]
            om = vorticity(base, points, tet, velocity)

            cur = current_liu_source(
                volume, q, d, om, nu, args.min_wall_distance, args.min_stilde
            )
            quad = quadrature_source(
                volume, q, d, om, nu, args.min_wall_distance, args.min_stilde
            )

            scale = (
                abs(cur["p_total"]) + abs(cur["d_total"])
                + abs(quad["p_total"]) + abs(quad["d_total"])
                + 1.0e-300
            )
            err = abs(
                (quad["p_total"] - quad["d_total"])
                - (cur["p_total"] - cur["d_total"])
            ) / scale

            for s in (total,):
                s["count"] += 1
                s["volume"] += volume
                s["p_current"] += cur["p_total"]
                s["d_current"] += cur["d_total"]
                s["p_quad"] += quad["p_total"]
                s["d_quad"] += quad["d_total"]
                s["err_sum"] += err
                s["err_max"] = max(s["err_max"], err)

            for s in bins:
                if s["d_lo"] <= dc < s["d_hi"]:
                    s["count"] += 1
                    s["volume"] += volume
                    s["p_current"] += cur["p_total"]
                    s["d_current"] += cur["d_total"]
                    s["p_quad"] += quad["p_total"]
                    s["d_quad"] += quad["d_total"]
                    s["err_sum"] += err
                    s["err_max"] = max(s["err_max"], err)
                    break

            if mu_e[ie] > 0.0:
                mut_ratio = mu_t_e[ie] / mu_e[ie]
                if mut_ratio > station_element_peak:
                    station_element_peak = mut_ratio
                    station_element_peak_y = yc

            worst.append((
                err, irank, ie + 1, xc, yc, dc, cur["qbar"] / nu,
                cur["p_total"], cur["d_total"],
                quad["p_total"], quad["d_total"],
                cur["r"], cur["fw"],
            ))

    xvals = sorted(set(
        round(v[0][0], 12)
        for v in owned_profile.values()
        if v[0][0] >= 0.0
    ))
    profile_x = min(xvals, key=lambda x: abs(x - args.station))
    station_nodes = [
        v for v in owned_profile.values()
        if abs(v[0][0] - profile_x) <= 2.0e-11
    ]
    nodal_peak = max((v[2] for v in station_nodes), default=float("nan"))
    nodal_peak_y = float("nan")
    for coord, _, ratio in station_nodes:
        if ratio == nodal_peak:
            nodal_peak_y = coord[1]
            break

    def serialise_stats(s):
        if s["count"] == 0:
            return {
                "d_range": [s["d_lo"], s["d_hi"]],
                "count": 0,
            }
        pc = s["p_current"]
        dc = s["d_current"]
        pq = s["p_quad"]
        dq = s["d_quad"]
        return {
            "d_range": [s["d_lo"], s["d_hi"]],
            "count": s["count"],
            "volume": s["volume"],
            "P_quad_over_current": safe_ratio(pq, pc),
            "D_quad_over_current": safe_ratio(dq, dc),
            "current_D_over_P": safe_ratio(dc, pc),
            "quad_D_over_P": safe_ratio(dq, pq),
            "current_net_P_minus_D": pc - dc,
            "quad_net_P_minus_D": pq - dq,
            "mean_normalised_net_difference": s["err_sum"] / s["count"],
            "max_normalised_net_difference": s["err_max"],
        }

    worst.sort(reverse=True)
    result = {
        "iteration": step,
        "station": args.station,
        "station_half_width": args.station_half_width,
        "station_profile_x": profile_x,
        "station_projected_nodal_mu_t_over_mu_peak": nodal_peak,
        "station_projected_nodal_peak_y": nodal_peak_y,
        "station_element_mu_t_over_mu_peak_in_window": station_element_peak,
        "station_element_peak_centroid_y": station_element_peak_y,
        "all_selected_elements": serialise_stats(total),
        "wall_distance_bins": [serialise_stats(s) for s in bins],
        "worst_elements": [
            {
                "normalised_net_difference": r[0],
                "rank": r[1],
                "local_element": r[2],
                "centroid_x": r[3],
                "centroid_y": r[4],
                "mean_wall_distance": r[5],
                "qbar_over_nu": r[6],
                "P_current": r[7],
                "D_current": r[8],
                "P_quad": r[9],
                "D_quad": r[10],
                "r_current": r[11],
                "fw_current": r[12],
            }
            for r in worst[:max(0, args.top)]
        ],
        "interpretation": {
            "D_quad_over_current_gt_1_outer":
                "current frozen-coefficient Liu source underestimates destruction in that band",
            "D_quad_over_current_lt_1_outer":
                "current frozen-coefficient Liu source overestimates destruction in that band",
            "near_unity":
                "four-point source quadrature is unlikely to explain the outer-layer error",
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=True))


if __name__ == "__main__":
    main()
