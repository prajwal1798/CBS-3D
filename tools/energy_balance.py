#!/usr/bin/env python3
# ================================
# FILE      : energy_balance.py
# PURPOSE   : STEADY-STATE MASS AND ENERGY BALANCE OF A CBS3D CHT RUN.
#             THE INLETS AND OUTLETS MUST LIE ON ONE PLANE (BLANKET SLICE:
#             y = 0). OVER EVERY TETRAHEDRON FACE ON THAT PLANE IT INTEGRATES
#               MASS      rho  * INT (u.n) dA
#               ENTHALPY  rho cp INT (u.n) (T - T_ref) dA   (EXACT FOR P1)
#             T_ref = MIXED-MEAN INLET TEMPERATURE, SO A SMALL MASS IMBALANCE
#             (rho cp dm T, LARGE AT T ~ 300 K) DOES NOT SWAMP THE HEAT.
#             WITH n THE OUTWARD NORMAL, SO INFLOW IS NEGATIVE. WETTED FACES
#             THAT SHARE A NODE FORM ONE PORT.
#             THE HEAT INPUT IS READ FROM THE SOLVER LOG LINE
#               HEAT_INPUT_TOTAL_W = ...
#             (OR GIVEN WITH --qin). CONVERGED WHEN Q_out / Q_in -> 1.
#
# USAGE     : python3 tools/energy_balance.py RESULT [--log solver.out]
#                 [--plane y=0] [--rho 997] [--cp 4182] [--qin W] [--tol 0.02]
#                 [--mean]
#             --mean USES THE TIME-AVERAGED FIELDS mean_velocity AND
#             mean_temperature OF A TRANSIENT RUN (CBS3D_AVERAGE_START_*).
#             WHEN THE FILE ALSO HAS covariance_velocity_temperature, IT
#             EVALUATES THE TRUE MEAN FLUX <u.n T> = <u.n><T> + <u.n'T'>.
#             WITHOUT THAT FIELD IT FALLS BACK TO <u.n><T>, WHICH OMITS THE
#             TURBULENT CORRELATION AND UNDERSTATES AN UNSTEADY OUTLET.
#             RESULT = .pvtu / .vtu / FOLDER (LATEST FILE UNDER IT)
#
# OUTPUT    : TABLE OF PORTS (MASS FLOW, MIXED-MEAN T), TOTALS, VERDICT.
#             EXIT 0 WHEN |Q_out/Q_in - 1| AND THE MASS IMBALANCE ARE BOTH
#             <= --tol (DEFAULT 0.02).
#
# STATUS    : NUMPY ONLY. CBS3D ASCII VTU.
# ================================

from __future__ import print_function

import argparse
import os
import re
import sys

import numpy as np

AXES = {"x": 0, "y": 1, "z": 2}
FACES = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])


# ================================
# CBS3D ASCII VTU READER
# ================================
def _array(text, name):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                  text, re.S)
    if m is None:
        raise ValueError("NO DataArray '%s'" % name)
    return np.array(m.group(1).split(), dtype=float)


def read_piece(path, velocity_name="velocity", temperature_name="temperature",
               covariance_name=None):
    with open(path) as handle:
        text = handle.read()
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    points = np.array(m.group(1).split(), dtype=float).reshape(-1, 3)
    conn = _array(text, "connectivity").astype(np.int64)
    offsets = _array(text, "offsets").astype(np.int64)
    types = _array(text, "types").astype(np.int64)
    tet = types == 10
    tets = conn[(offsets[tet] - 4)[:, None] + np.arange(4)[None, :]]
    velocity = _array(text, velocity_name).reshape(-1, 3)
    temperature = _array(text, temperature_name)
    covariance = None
    if covariance_name is not None:
        try:
            covariance = _array(text, covariance_name).reshape(-1, 3)
        except ValueError:
            covariance = None       # OLDER OUTPUT: PRODUCT <u><T> ONLY
    return points, tets, velocity, temperature, covariance


def latest(path):
    if not os.path.isdir(path):
        return path
    found = []
    for root, _dirs, files in os.walk(path):
        for f in files:
            if f.endswith(".pvtu") or (f.endswith(".vtu") and "_rank_" not in f
                                       and "boundary" not in f):
                full = os.path.join(root, f)
                found.append((os.path.getmtime(full), full))
    if not found:
        raise ValueError("NO .pvtu OR .vtu UNDER %s" % path)
    return max(found)[1]


def read_result(path, mean=False):
    path = latest(path)
    if mean:
        names = ("mean_velocity", "mean_temperature")
        covariance_name = "covariance_velocity_temperature"
    else:
        names = ("velocity", "temperature")
        covariance_name = None
    if path.endswith(".pvtu"):
        with open(path) as handle:
            sources = re.findall(r'<Piece\s+Source="([^"]+)"', handle.read())
        base = os.path.dirname(path)
        return [read_piece(os.path.join(base, s), *names,
                           covariance_name=covariance_name) for s in sources], path
    return [read_piece(path, *names, covariance_name=covariance_name)], path


def heat_input_from_log(log):
    value = None
    with open(log) as handle:
        for line in handle:
            m = re.search(r"HEAT_INPUT_TOTAL_W\s*=\s*([-+0-9.eE]+)", line)
            if m:
                value = float(m.group(1))
    return value


# ================================
# FACES ON THE PORT PLANE
# ================================
def plane_faces(pieces, axis, value, tol):
    """AREA, u.e_axis, T, cov(u_axis,T), CENTROID, NODE KEYS."""
    areas, un, temp, cov, cent, keys = [], [], [], [], [], []
    use_cov = all(piece[4] is not None for piece in pieces)
    for ip, (points, tets, vel, t, piece_cov) in enumerate(pieces):
        on = np.abs(points[:, axis] - value) <= tol
        for f in FACES:
            tri = tets[:, f]
            keep = on[tri].all(axis=1)
            if not keep.any():
                continue
            tri = tri[keep]
            p0, p1, p2 = points[tri[:, 0]], points[tri[:, 1]], points[tri[:, 2]]
            areas.append(0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1))
            un.append(vel[tri][:, :, axis])
            temp.append(t[tri])
            if use_cov:
                cov.append(piece_cov[tri][:, :, axis])
            cent.append((p0 + p1 + p2) / 3.0)
            # NODES ARE MATCHED ACROSS PIECES BY THEIR ROUNDED COORDINATES
            keys.append(np.round(points[tri], 9).reshape(-1, 9))
    if not areas:
        raise ValueError("NO TETRAHEDRON FACE LIES ON THE PLANE")
    return (np.concatenate(areas), np.vstack(un), np.vstack(temp),
            np.vstack(cov) if use_cov else None,
            np.vstack(cent), np.vstack(keys))


def group_ports(keys):
    """PORT INDEX PER FACE: FACES SHARING A NODE BELONG TO ONE PORT."""
    parent = list(range(len(keys)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    owner = {}
    for f, row in enumerate(keys):
        for n in range(3):
            node = tuple(row[3 * n:3 * n + 3])
            if node in owner:
                a, b = find(f), find(owner[node])
                if a != b:
                    parent[a] = b
            else:
                owner[node] = f
    roots = {}
    return np.array([roots.setdefault(find(f), len(roots)) for f in range(len(keys))])


def main():
    ap = argparse.ArgumentParser(description="CBS3D mass and energy balance")
    ap.add_argument("result")
    ap.add_argument("--log", help="SOLVER LOG WITH HEAT_INPUT_TOTAL_W")
    ap.add_argument("--qin", type=float, help="HEAT INPUT [W] (OVERRIDES --log)")
    ap.add_argument("--plane", default="y=0", help="PORT PLANE, E.G. y=0")
    ap.add_argument("--rho", type=float, default=997.0)
    ap.add_argument("--cp", type=float, default=4182.0)
    ap.add_argument("--tref", type=float,
                    help="REFERENCE T [K] (DEFAULT: MIXED-MEAN INLET T)")
    ap.add_argument("--tol", type=float, default=0.02,
                    help="PASS WHEN |Q_out/Q_in - 1| <= tol")
    ap.add_argument("--mean", action="store_true",
                    help="USE THE TIME-AVERAGED FIELDS (mean_velocity, mean_temperature)")
    args = ap.parse_args()

    name, val = args.plane.split("=")
    axis = AXES[name.strip().lower()]
    value = float(val)

    try:
        pieces, path = read_result(args.result, mean=args.mean)
    except ValueError as error:
        if args.mean:
            print("ERROR   : %s; THIS FILE HAS NO TIME AVERAGES (CBS3D_AVERAGE_START_* "
                  "NOT SET, OR THE WINDOW HAD NOT STARTED)" % error)
            return 2
        raise
    print("READ    : %s (%d piece(s))%s" % (path, len(pieces),
                                         ", TIME-AVERAGED FIELDS" if args.mean else ""))

    allpts = np.vstack([p[0] for p in pieces])
    extent = allpts.max(axis=0) - allpts.min(axis=0)
    tol = 1.0e-7 * max(extent.max(), 1.0e-12)
    lo, hi = allpts[:, axis].min(), allpts[:, axis].max()
    sign = -1.0 if abs(value - lo) <= abs(value - hi) else 1.0   # OUTWARD

    area, un, temp, cov, cent, keys = plane_faces(pieces, axis, value, tol)
    un = sign * un                                                  # u.n_out
    if cov is not None:
        cov = sign * cov                                            # cov(u.n,T)
    def enthalpy(t_shift):
        tt = temp - t_shift
        return args.rho * args.cp * area / 12.0 * (
            un.sum(axis=1) * tt.sum(axis=1) + (un * tt).sum(axis=1))

    def enthalpy_with_covariance(t_shift):
        # INT u.n (T - T_ref) = A/12[(sum u)(sum T) + sum <uT>]
        #                      - T_ref A/3 sum u,   <uT> = <u><T> + cov
        nodal_ut = un * temp + cov
        return args.rho * args.cp * (
            area / 12.0 * (un.sum(axis=1) * temp.sum(axis=1) +
                           nodal_ut.sum(axis=1)) -
            t_shift * area / 3.0 * un.sum(axis=1))

    flux = enthalpy_with_covariance if cov is not None else enthalpy

    mass = args.rho * area / 3.0 * un.sum(axis=1)
    inflow = mass < 0.0
    t_in = flux(0.0)[inflow].sum() / (args.cp * mass[inflow].sum()) \
        if inflow.any() else 0.0
    t_ref = args.tref if args.tref is not None else t_in
    enth = flux(t_ref)
    q_cov = enth.sum() - enthalpy(t_ref).sum() if cov is not None else 0.0

    if args.mean:
        if cov is None:
            print("FLUX    : rho cp INT <u.n><T> dA; covariance_velocity_temperature "
                  "IS ABSENT")
        else:
            print("FLUX    : rho cp INT <u.n T> dA = <u.n><T> + <u.n'T'> "
                  "(covariance field present)")

    # PORTS: CONNECTED GROUPS OF WETTED FACES, LISTED ALONG THE IN-PLANE
    # DIRECTION OF LARGEST EXTENT
    inplane = [a for a in range(3) if a != axis]
    along = max(inplane, key=lambda a: extent[a])
    wet = np.abs(un).max(axis=1) > 0.0
    port = group_ports(keys[wet])
    order = np.argsort([np.mean(cent[wet, along][port == k])
                        for k in range(port.max() + 1)])
    port = np.argsort(order)[port]

    print("PLANE   : %s = %g, outward normal %+d, %d faces (%d wetted)"
          % ("xyz"[axis], value, int(sign), area.size, int(wet.sum())))
    print()
    print("T_REF   : %.4f K (%s)" % (t_ref, "--tref" if args.tref is not None
                                       else "mixed-mean inlet"))
    print("%-5s %12s %14s %14s %12s" % ("PORT", "%s_mean [m]" % "xyz"[along],
                                       "mdot [kg/s]", "Q [W]", "T_mix [K]"))
    m_w, e_w, c_w, a_w = mass[wet], enth[wet], cent[wet, along], area[wet]
    for k in range(port.max() + 1):
        sel = port == k
        mk, ek = m_w[sel].sum(), e_w[sel].sum()
        tmix = t_ref + ek / (mk * args.cp) if mk != 0.0 else float("nan")
        kind = "out" if mk > 0 else "in"
        print("%-5s %12.5f %14.6e %14.6e %12.4f"
              % ("%d%s" % (k + 1, kind[0]), np.average(c_w[sel], weights=a_w[sel]),
                 mk, ek, tmix))

    m_in = -mass[mass < 0].sum()
    m_out = mass[mass > 0].sum()
    q_out = enth.sum()
    print()
    print("MASS    : in %.6e  out %.6e kg/s  (imbalance %.3e of inflow)"
          % (m_in, m_out, (m_out - m_in) / max(m_in, 1e-300)))
    print("HEAT    : carried away by the coolant Q_out = %.6e W"
          "  (sum of mdot cp (T - T_ref))" % q_out)
    if args.mean and cov is not None:
        print("COVARIANCE: <u.n'T'> term = %+.6e W (%.2f%% of Q_out)"
              % (q_cov, 100.0 * q_cov / q_out if q_out != 0.0 else float("nan")))
    elif args.mean:
        print("NOTE    : THE PRODUCT <u.n><T> OMITS <u.n'T'>; AN UNSTEADY OUTLET "
              "CAN SHOW A SPURIOUS DEFICIT")
    t_all = np.concatenate([p[3] for p in pieces])
    print("T RANGE : %.4f .. %.4f K (all nodes)" % (t_all.min(), t_all.max()))

    q_in = args.qin
    if q_in is None and args.log:
        q_in = heat_input_from_log(args.log)
    if q_in is None:
        print("Q_in    : UNKNOWN (GIVE --log solver.out OR --qin)")
        return 2
    ratio = q_out / q_in if q_in != 0.0 else float("nan")
    print("Q_in    : %.6e W   Q_out / Q_in = %.5f" % (q_in, ratio))
    m_err = abs(m_out - m_in) / max(m_in, 1e-300)
    e_ok = abs(ratio - 1.0) <= args.tol
    m_ok = m_err <= args.tol
    print("ENERGY  : %s (|Q_out/Q_in - 1| = %.3g, limit %.3g)"
          % ("BALANCED" if e_ok else "NOT BALANCED YET", abs(ratio - 1.0), args.tol))
    print("MASS    : %s (|imbalance| = %.3g, limit %.3g)%s"
          % ("BALANCED" if m_ok else "NOT BALANCED", m_err, args.tol,
             "" if m_ok else "; IF IT STAYS CONSTANT IT IS THE DISCRETE "
             "CONTINUITY ERROR OF THE MESH, NOT A TRANSIENT"))
    if e_ok and m_ok:
        print("RESULT  : STEADY")
        return 0
    print("RESULT  : NOT STEADY / NOT BALANCED")
    if args.mean and cov is None:
        print("CAVEAT  : THIS MEAN FIELD HAS NO <u'T'> DATA; THE DEFICIT ABOVE "
              "IS NOT A SOLVER BALANCE")
    return 1


if __name__ == "__main__":
    sys.exit(main())
