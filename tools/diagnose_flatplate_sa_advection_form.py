#!/usr/bin/env python3
"""Frozen-field audit of conservative versus non-conservative SA advection.

The incompressible SA equation is normally written with u.grad(nu_tilde).
The current Liu production assembly uses div(u nu_tilde), which differs by
nu_tilde div(u).  These forms are identical only for pointwise divergence-free
velocity.  P1 CBS velocity satisfies incompressibility only discretely, so this
script measures the extra term on an existing distributed flat-plate snapshot.

No CFD solve is performed and no third-party Python package is required.
"""
from __future__ import print_function

import argparse
import importlib.util
import json
import math
from pathlib import Path

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


def load_base(repo_root):
    path = repo_root / "examples" / "Flat_Plate_Turbulent" / "analyze_distributed_tmr.py"
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
    p.add_argument("--top", type=int, default=12)
    return p.parse_args()


def tet_volume(base, points, tet):
    p0, p1, p2, p3 = (points[n - 1] for n in tet)
    return abs(base.dot(base.sub(p1, p0), base.cross(base.sub(p2, p0), base.sub(p3, p0)))) / 6.0


def new_stats(lo, hi):
    return {
        "d_lo": lo, "d_hi": hi, "count": 0, "volume": 0.0,
        "abs_nonconservative": 0.0, "abs_conservative": 0.0,
        "abs_divergence_correction": 0.0,
        "signed_divergence_correction": 0.0,
        "abs_div_u_volume": 0.0, "signed_div_u_volume": 0.0,
        "max_abs_div_u": 0.0, "max_normalised_div_u": 0.0,
        "normalised_div_u_sum": 0.0,
    }


def serialise(s):
    if s["count"] == 0:
        return {"d_range": [s["d_lo"], s["d_hi"]], "count": 0}
    eps = 1.0e-300
    return {
        "d_range": [s["d_lo"], s["d_hi"]],
        "count": s["count"],
        "volume": s["volume"],
        "abs_divergence_correction_over_abs_nonconservative":
            s["abs_divergence_correction"] / max(s["abs_nonconservative"], eps),
        "abs_conservative_over_abs_nonconservative":
            s["abs_conservative"] / max(s["abs_nonconservative"], eps),
        "signed_divergence_correction_over_abs_nonconservative":
            s["signed_divergence_correction"] / max(s["abs_nonconservative"], eps),
        "volume_mean_abs_div_u": s["abs_div_u_volume"] / max(s["volume"], eps),
        "volume_mean_signed_div_u": s["signed_div_u_volume"] / max(s["volume"], eps),
        "max_abs_div_u": s["max_abs_div_u"],
        "mean_abs_div_u_over_grad_u_frobenius":
            s["normalised_div_u_sum"] / float(s["count"]),
        "max_abs_div_u_over_grad_u_frobenius": s["max_normalised_div_u"],
    }


def main():
    args = parse_args()
    base = load_base(args.repo_root)
    ranks = sorted(p for p in args.partition_root.glob("rank_[0-9][0-9][0-9][0-9]") if p.is_dir())
    if len(ranks) != 40:
        raise SystemExit("expected 40 partition ranks, found {}".format(len(ranks)))

    rx = base.re.compile(r"flatplate_step_(\d{8})_rank_(\d{4})\.vtu$")
    by_step = {}
    for f in args.output_dir.glob("flatplate_step_*_rank_*.vtu"):
        m = rx.match(f.name)
        if m:
            by_step.setdefault(int(m.group(1)), {})[int(m.group(2))] = f
    good = [step for step, pieces in by_step.items() if len(pieces) == 40]
    if not good:
        raise SystemExit("no complete 40-piece distributed snapshot")
    step = max(good)

    bins = [new_stats(lo, hi) for lo, hi in BINS]
    total = new_stats(0.0, float("inf"))
    worst = []

    for irank, rank in enumerate(ranks):
        tets, mesh_points, _, _ = base.read_partition(rank)
        points, velocity, pd, cd = base.read_piece(by_step[step][irank])
        qfield = pd["nu_tilde"]
        dfield = pd["wall_distance"]
        mu_e = cd.get("mu_e", [])

        for ie, tet in enumerate(tets):
            if mu_e and not (mu_e[ie] > 0.0):
                continue
            xc = 0.25 * sum(points[n - 1][0] for n in tet)
            if abs(xc - args.station) > args.station_half_width:
                continue
            dc = 0.25 * sum(dfield[n - 1] for n in tet)
            volume = tet_volume(base, points, tet)
            grads = base.tet_gradients(points, tet)

            q = [qfield[n - 1] for n in tet]
            grad_q = [0.0, 0.0, 0.0]
            for a, node in enumerate(tet):
                for j in range(3):
                    grad_q[j] += q[a] * grads[a][j]

            gu = base.velocity_gradient(points, tet, velocity)
            div_u = gu[0][0] + gu[1][1] + gu[2][2]
            grad_u_frob = math.sqrt(sum(gu[i][j] * gu[i][j] for i in range(3) for j in range(3)))
            normalised_div = abs(div_u) / max(grad_u_frob, 1.0e-300)

            usum = [0.0, 0.0, 0.0]
            for node in tet:
                v = velocity[node - 1]
                for j in range(3):
                    usum[j] += v[j]

            nc = [0.0] * 4
            corr = [0.0] * 4
            cons = [0.0] * 4
            qsum = sum(q)
            for a, node in enumerate(tet):
                va = velocity[node - 1]
                weighted_u_dot_gradq = sum((usum[j] + va[j]) * grad_q[j] for j in range(3))
                nc[a] = -(volume / 20.0) * weighted_u_dot_gradq
                qload = volume * (qsum + q[a]) / 20.0
                corr[a] = -div_u * qload
                cons[a] = nc[a] + corr[a]

            abs_nc = sum(abs(v) for v in nc)
            abs_cons = sum(abs(v) for v in cons)
            abs_corr = sum(abs(v) for v in corr)
            signed_corr = sum(corr)

            targets = [total]
            for s in bins:
                if s["d_lo"] <= dc < s["d_hi"]:
                    targets.append(s)
                    break
            for s in targets:
                s["count"] += 1
                s["volume"] += volume
                s["abs_nonconservative"] += abs_nc
                s["abs_conservative"] += abs_cons
                s["abs_divergence_correction"] += abs_corr
                s["signed_divergence_correction"] += signed_corr
                s["abs_div_u_volume"] += abs(div_u) * volume
                s["signed_div_u_volume"] += div_u * volume
                s["max_abs_div_u"] = max(s["max_abs_div_u"], abs(div_u))
                s["normalised_div_u_sum"] += normalised_div
                s["max_normalised_div_u"] = max(s["max_normalised_div_u"], normalised_div)

            ratio = abs_corr / max(abs_nc, 1.0e-300)
            yc = 0.25 * sum(points[n - 1][1] for n in tet)
            worst.append((ratio, irank, ie + 1, xc, yc, dc, div_u, normalised_div, abs_nc, abs_corr))

    worst.sort(reverse=True)
    result = {
        "iteration": step,
        "station": args.station,
        "station_half_width": args.station_half_width,
        "meaning": {
            "divergence_correction": "current conservative SA advection minus incompressible non-conservative SA advection",
            "ratio_near_zero": "conservative and non-conservative forms are numerically equivalent on this field",
            "ratio_order_one": "q div(u) materially changes SA advection and must be treated as a leading suspect",
        },
        "all_selected_elements": serialise(total),
        "wall_distance_bins": [serialise(s) for s in bins],
        "worst_elements": [
            {
                "abs_correction_over_abs_nonconservative": r[0],
                "rank": r[1], "local_element": r[2],
                "centroid_x": r[3], "centroid_y": r[4], "mean_wall_distance": r[5],
                "div_u": r[6], "abs_div_u_over_grad_u_frobenius": r[7],
                "abs_nonconservative_advection": r[8], "abs_divergence_correction": r[9],
            }
            for r in worst[:max(1, args.top)]
        ],
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
