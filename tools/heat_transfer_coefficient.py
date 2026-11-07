#!/usr/bin/env python3
# ================================
# FILE      : heat_transfer_coefficient.py
# PURPOSE   : WALL HEAT FLUX, HEAT TRANSFER COEFFICIENT AND NUSSELT NUMBER ON
#             THE COOLANT-STEEL INTERFACE OF A CBS3D CONJUGATE RESULT.
#
#   INTERFACE   EVERY TETRAHEDRON FACE SHARED BY A FLUID (material 0) AND A
#               SOLID ELEMENT. n = UNIT NORMAL FROM THE SOLID INTO THE FLUID.
#   WALL FLUX   CONSISTENT (RESIDUAL) FLUX OF THE STEEL: AT EVERY INTERFACE
#               NODE a THE HEAT LEAVING THE STEEL IS
#                 F_a = sum_solid [ V Q N_a  +  BC532 LOAD  -  k V grad N_a . grad T ]
#               (THE SOLID-ONLY FINITE-ELEMENT EQUATION OF THAT NODE), SPREAD
#               OVER ITS LUMPED INTERFACE AREA: q_w = F_a / A_a. THE INTEGRAL
#               OF q_w THEN EQUALS THE HEAT INPUT EXACTLY AT STEADY STATE,
#               WHICH ELEMENT GRADIENTS DO NOT GUARANTEE (FLUID-SIDE P1
#               GRADIENTS MISS MOST OF A THIN THERMAL BOUNDARY LAYER).
#               THE ELEMENT-GRADIENT VALUES -k grad(T) . n ARE ALSO WRITTEN
#               (q_w_solid_gradient, q_w_fluid_gradient). POSITIVE = INTO THE
#               COOLANT. FOR A TRANSIENT SNAPSHOT THE CONSISTENT FLUX OMITS THE
#               STEEL'S HEAT STORAGE; USE IT FOR (QUASI-)STEADY STATES.
#   HEATING     FROM --par (THE RUN'S .par: heat_flux_bc, source_solid AND
#               THE NUCLEAR-HEATING PROFILES) OR, WITHOUT IT, THE BLANKET L2
#               VALUES; --heat-scale MULTIPLIES BOTH (PULSE FACTOR f(t)).
#   WALL T      T_w = MEAN OF THE THREE FACE-NODE TEMPERATURES
#   h_in        q_w / (T_w - T_in)           T_in = INLET TEMPERATURE
#   h_b         q_w / (T_w - T_b)            T_b  = LOCAL BULK (MIXING-CUP)
#               TEMPERATURE: VOLUME- AND u.e-WEIGHTED MEAN OVER THE COOLANT OF
#               THE SAME CHANNEL IN A THIN SLAB |(x - x_f).e| < 0.15 D_h
#               NORMAL TO THE LOCAL MEAN FLOW DIRECTION e, WITHIN 1.6 D_h OF
#               THE FACE. APPROXIMATE IN BENDS, WHERE NO CROSS-SECTION IS PLANE.
#   Nu          h D_h / k_f   (Nu_in, Nu_b)
#
#   CHANNELS    CONNECTED COOLANT REGIONS, NUMBERED BY DECREASING MAXIMUM y
#               (BLANKET SLICE: 1 = NEXT TO THE FIRST WALL). --dh GIVES THEIR
#               HYDRAULIC DIAMETERS IN THAT ORDER (SQUARE CHANNELS: D_h = WIDTH).
#
#   CHECK       AT STEADY STATE ALL HEAT LEAVES THROUGH THE COOLANT, SO THE
#               INTEGRAL OF q_w OVER THE INTERFACE MUST EQUAL THE HEAT INPUT
#               (HEAT_INPUT_TOTAL_W IN THE LOG, TIMES THE PULSE FACTOR).
#
# OUTPUT    : <out>_interface.vtu   TRIANGLES WITH q_w, q_w_fluid, T_w, T_b,
#                                   h_in, h_b, Nu_in, Nu_b, channel (PARAVIEW)
#             <out>_channels.csv    PER CHANNEL: AREA, HEAT, MEAN h, MEAN Nu
#             <out>_profiles.png    q_w, T_w AND T_b, h_b, Nu_b ALONG THE
#                                   HEATED (FIRST-WALL-SIDE) WALL OF EACH CHANNEL
#             <out>_profiles.csv    THE SAME DATA
#
# USAGE     : python3 tools/heat_transfer_coefficient.py RESULT
#                 [--out PREFIX] [--k-fluid 0.606] [--k-solid 28.3]
#                 [--t-in 298.15] [--dh 0.015 0.010 0.010 0.005]
#                 [--log solver.out | --qin W] [--axis x]
#                 [--par RUN.par] [--heat-scale f]
#             RESULT = .pvtu / .vtu / FOLDER (LATEST FILE UNDER IT)
#
# STATUS    : NUMPY + SCIPY (cKDTree, connected components) + MATPLOTLIB.
#             CBS3D ASCII VTU (global_node_id, material_id).
# ================================

from __future__ import print_function

import argparse
import csv
import os
import re
import sys

import numpy as np

FACES = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])   # FACE OPPOSITE NODE i


# ================================
# READER (ALL PIECES -> ONE GLOBAL MESH BY global_node_id)
# ================================
def _array(text, name):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, text, re.S)
    if m is None:
        raise ValueError("NO DataArray '%s'" % name)
    return np.array(m.group(1).split(), dtype=float)


def latest(path):
    if not os.path.isdir(path):
        return path
    found = []
    for root, _d, files in os.walk(path):
        for f in files:
            if f.endswith(".pvtu"):
                full = os.path.join(root, f)
                found.append((os.path.getmtime(full), full))
    if not found:
        raise ValueError("NO .pvtu UNDER %s" % path)
    return max(found)[1]


def read_mesh(path, velocity_name, temperature_name="temperature"):
    path = latest(path)
    if path.endswith(".pvtu"):
        with open(path) as handle:
            sources = re.findall(r'<Piece\s+Source="([^"]+)"', handle.read())
        files = [os.path.join(os.path.dirname(path), s) for s in sources]
    else:
        files = [path]
    xs, ts, us, gids, tets, mats, tris, bcs = [], [], [], [], [], [], [], []
    for f in files:
        with open(f) as handle:
            text = handle.read()
        m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
        x = np.array(m.group(1).split(), dtype=float).reshape(-1, 3)
        gid = _array(text, "global_node_id").astype(np.int64)
        conn = _array(text, "connectivity").astype(np.int64)
        offsets = _array(text, "offsets").astype(np.int64)
        types = _array(text, "types").astype(np.int64)
        mat = _array(text, "material_id").astype(np.int64)
        tet = types == 10
        local = conn[(offsets[tet] - 4)[:, None] + np.arange(4)[None, :]]
        tets.append(gid[local])
        mats.append(mat[tet])
        tri = types == 5
        if tri.any():
            tris.append(gid[conn[(offsets[tri] - 3)[:, None] + np.arange(3)[None, :]]])
            bcs.append(_array(text, "bc_id").astype(np.int64)[tri])
        xs.append(x)
        gids.append(gid)
        ts.append(_array(text, temperature_name))
        us.append(_array(text, velocity_name).reshape(-1, 3))
    gid = np.concatenate(gids)
    n = gid.max() + 1
    X = np.zeros((n, 3))
    T = np.zeros(n)
    U = np.zeros((n, 3))
    X[gid] = np.concatenate(xs)
    T[gid] = np.concatenate(ts)
    U[gid] = np.concatenate(us)
    used = np.zeros(n, bool)
    used[gid] = True
    btri = np.concatenate(tris) if tris else np.zeros((0, 3), np.int64)
    bc = np.concatenate(bcs) if bcs else np.zeros(0, np.int64)
    return path, X, T, U, used, np.concatenate(tets), np.concatenate(mats), btri, bc


# ================================
# HEATING FROM A CBS3D .par (OR THE BLANKET L2 DEFAULTS)
# ================================
def read_heating(par):
    heat = dict(flux=1.0e5, source_solid=0.0, nuclear_on=1, origin=(0.0, 0.28, 0.0),
                normal=(0.0, -1.0, 0.0), profiles=[(1, 1.6e6, 0.15)], source="blanket L2 defaults")
    if par is None:
        return heat
    lines = [l.rstrip("\n") for l in open(par) if not l.lstrip().startswith("#")]
    for i, l in enumerate(lines):
        if l.startswith("alpha_sf k_ratio source_solid heat_flux_bc"):
            v = lines[i + 1].split()
            heat["source_solid"], heat["flux"] = float(v[2]), float(v[3])
        if l.startswith("Nuclear heating: nuclear_heating_on"):
            v = lines[i + 1].split()
            heat["nuclear_on"] = int(v[0])
            heat["origin"] = tuple(map(float, v[1:4]))
            heat["normal"] = tuple(map(float, v[4:7]))
            nprof = int(v[7])
            heat["profiles"] = [(int(lines[i + 3 + k].split()[0]), float(lines[i + 3 + k].split()[1]),
                                 float(lines[i + 3 + k].split()[2])) for k in range(nprof)]
    if not any(l.startswith("Nuclear heating: nuclear_heating_on") for l in lines):
        heat["nuclear_on"] = 0
    heat["source"] = par
    return heat


def find_par(result_path):
    """<run>/partitions/rank_0000/*_rank_0000.par FOR A RESULT IN <run>/output/<case>/."""
    run = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(result_path))))
    folder = os.path.join(run, "partitions", "rank_0000")
    if os.path.isdir(folder):
        for f in sorted(os.listdir(folder)):
            if f.endswith(".par"):
                return os.path.join(folder, f)
    return None


def heat_input_from_log(log):
    value, scale = None, 1.0
    with open(log) as handle:
        for line in handle:
            m = re.search(r"HEAT_INPUT_TOTAL_W\s*=\s*([-+0-9.eE]+)", line)
            if m:
                value = float(m.group(1))
    return value


# ================================
# P1 GRADIENT OF T IN A SET OF TETRAHEDRA
# ================================
def tet_gradient(X, T, tets):
    x0 = X[tets[:, 0]]
    A = np.stack([X[tets[:, 1]] - x0, X[tets[:, 2]] - x0, X[tets[:, 3]] - x0], axis=1)
    b = np.stack([T[tets[:, 1]] - T[tets[:, 0]], T[tets[:, 2]] - T[tets[:, 0]],
                  T[tets[:, 3]] - T[tets[:, 0]]], axis=1)
    return np.linalg.solve(A, b[:, :, None])[:, :, 0]


def tet_volume(X, tets):
    x0 = X[tets[:, 0]]
    return np.abs(np.einsum("ij,ij->i", X[tets[:, 1]] - x0,
                            np.cross(X[tets[:, 2]] - x0, X[tets[:, 3]] - x0))) / 6.0


# ================================
# MAIN
# ================================
def main():
    ap = argparse.ArgumentParser(description=__doc__ if __doc__ else "")
    ap.add_argument("result")
    ap.add_argument("--out", default=None)
    ap.add_argument("--velocity", default="velocity", help="velocity array (e.g. mean_velocity)")
    ap.add_argument("--k-fluid", type=float, default=0.606)
    ap.add_argument("--k-solid", type=float, default=28.3)
    ap.add_argument("--t-in", type=float, default=298.15)
    ap.add_argument("--dh", type=float, nargs="+", default=[0.015, 0.010, 0.010, 0.005])
    ap.add_argument("--log", default=None)
    ap.add_argument("--qin", type=float, default=None, help="heat input [W] (instead of --log)")
    ap.add_argument("--axis", default="x", help="axis for the profile plots along the heated walls")
    ap.add_argument("--par", default=None, help="the run's .par (heating); default: found next to the run")
    ap.add_argument("--heat-scale", type=float, default=1.0, help="pulse factor f(t) of this snapshot")
    ap.add_argument("--bin", type=float, default=1.0, help="profile bin width in units of D_h")
    args = ap.parse_args()

    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    path, X, T, U, used, tets, mats, btri, bc = read_mesh(args.result, args.velocity)
    heat = read_heating(args.par or find_par(path))
    print("HEATING: %s: q'' = %g W/m2, nuclear %s, scale f = %g"
          % (heat["source"], heat["flux"], heat["profiles"] if heat["nuclear_on"] else "off", args.heat_scale))
    out = args.out or os.path.splitext(path)[0] + "_htc"
    print("RESULT : %s (%d tetrahedra)" % (path, len(tets)))

    # ================================
    # FACES: SORT NODE TRIPLETS, PAIR EQUAL ONES
    # ================================
    nt = len(tets)
    faces = np.sort(tets[:, FACES].reshape(-1, 3), axis=1)            # (4 nt, 3)
    owner = np.repeat(np.arange(nt), 4)
    opposite = tets[:, [0, 1, 2, 3]].reshape(-1)                      # NODE OPPOSITE EACH FACE
    order = np.lexsort((faces[:, 2], faces[:, 1], faces[:, 0]))
    fs = faces[order]
    same = np.all(fs[1:] == fs[:-1], axis=1)
    i1 = order[:-1][same]
    i2 = order[1:][same]                                              # INTERIOR FACE PAIRS
    t1, t2 = owner[i1], owner[i2]
    fluid = mats == 0

    # COOLANT CHANNELS = CONNECTED FLUID REGIONS
    both = fluid[t1] & fluid[t2]
    adj = coo_matrix((np.ones(both.sum()), (t1[both], t2[both])), shape=(nt, nt))
    ncomp, comp = connected_components(adj, directed=False)
    fcomp = np.unique(comp[fluid])
    ymax = np.array([X[tets[comp == c]].reshape(-1, 3)[:, 1].max() for c in fcomp])
    rank = fcomp[np.argsort(-ymax)]
    channel_of_comp = {c: k + 1 for k, c in enumerate(rank)}
    nch = len(rank)
    if len(args.dh) < nch:
        sys.exit("ERROR: %d coolant channels found but only %d values in --dh" % (nch, len(args.dh)))
    print("CHANNELS: %d (by decreasing max y), D_h = %s m" % (nch, " ".join("%g" % d for d in args.dh[:nch])))

    # ================================
    # INTERFACE FACES: ONE FLUID AND ONE SOLID TET
    # ================================
    mixed = fluid[t1] != fluid[t2]
    fi = np.where(fluid[t1[mixed]], i1[mixed], i2[mixed])            # FACE INDEX ON THE FLUID SIDE
    si = np.where(fluid[t1[mixed]], i2[mixed], i1[mixed])            # ... ON THE SOLID SIDE
    tf, ts_ = owner[fi], owner[si]
    tri = faces[si]
    nf = len(tri)
    a, b, c = X[tri[:, 0]], X[tri[:, 1]], X[tri[:, 2]]
    n = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(n, axis=1)
    n /= (2.0 * area)[:, None]
    # ORIENT n OUT OF THE SOLID TET (AWAY FROM ITS OPPOSITE NODE) = INTO THE FLUID
    d = X[opposite[si]] - a
    flip = np.einsum("ij,ij->i", n, d) > 0
    n[flip] *= -1.0
    centre = (a + b + c) / 3.0
    gs = tet_gradient(X, T, tets[ts_])
    gf = tet_gradient(X, T, tets[tf])
    q_sg = -args.k_solid * np.einsum("ij,ij->i", gs, n)
    q_f = -args.k_fluid * np.einsum("ij,ij->i", gf, n)

    # ================================
    # CONSISTENT FLUX FROM THE STEEL'S NODAL EQUATIONS
    # ================================
    st = tets[~fluid]
    M = np.concatenate([np.ones((len(st), 4, 1)), X[st]], axis=2)       # [1 x y z] PER NODE
    gradN = np.linalg.inv(M)[:, 1:, :]                                  # (ns, 3, 4): dN_a/dx_j
    Vs = tet_volume(X, st)
    gT = np.einsum("eja,ea->ej", gradN, T[st])
    conduction = args.k_solid * Vs[:, None] * np.einsum("eja,ej->ea", gradN, gT)
    Qe = np.full(len(st), heat["source_solid"])
    if heat["nuclear_on"]:
        nn_ = np.array(heat["normal"]) / np.linalg.norm(heat["normal"])
        depth = np.maximum(0.0, (X[st].mean(axis=1) - np.array(heat["origin"])) @ nn_)
        smat = mats[~fluid]
        for mid, A_, lam in heat["profiles"]:
            Qe = np.where(smat == mid, A_ * np.exp(-depth / lam), Qe)
    F = np.zeros(len(X))
    for k in range(4):
        np.add.at(F, st[:, k], args.heat_scale * Qe * Vs / 4.0 - conduction[:, k])
    heated = btri[bc == 532]
    if len(heated):
        ha = 0.5 * np.linalg.norm(np.cross(X[heated[:, 1]] - X[heated[:, 0]],
                                           X[heated[:, 2]] - X[heated[:, 0]]), axis=1)
        for k in range(3):
            np.add.at(F, heated[:, k], args.heat_scale * heat["flux"] * ha / 3.0)
    on_interface = np.zeros(len(X), bool)
    on_interface[tri.reshape(-1)] = True
    A_node = np.zeros(len(X))
    for k in range(3):
        np.add.at(A_node, tri[:, k], area / 3.0)
    q_node = np.where(A_node > 0, F / np.where(A_node > 0, A_node, 1.0), 0.0)
    q_s = q_node[tri].mean(axis=1)
    off = ~on_interface & used
    solid_nodes = np.zeros(len(X), bool)
    solid_nodes[st.reshape(-1)] = True
    resid = np.abs(F[off & solid_nodes]).sum()
    T_w = T[tri].mean(axis=1)
    ch = np.array([channel_of_comp[c_] for c_ in comp[tf]])
    dh = np.array(args.dh)[ch - 1]
    print("INTERFACE: %d faces, area %.6g m2" % (nf, area.sum()))

    # ================================
    # LOCAL BULK TEMPERATURE PER CHANNEL (SLAB MIXING CUP)
    # ================================
    vol = tet_volume(X, tets)
    node_vol = np.zeros(len(X))
    for k in range(4):
        np.add.at(node_vol, tets[fluid, k], vol[fluid] / 4.0)
    T_b = np.full(nf, np.nan)
    for k in range(1, nch + 1):
        comp_id = rank[k - 1]
        nodes = np.unique(tets[(comp == comp_id)].reshape(-1))
        tree = cKDTree(X[nodes])
        sel = np.where(ch == k)[0]
        D = args.dh[k - 1]
        # ONE STATION PER 0.25 D_h CUBE OF FACE CENTRES; FACES SHARE ITS T_b
        key = np.floor(centre[sel] / (0.25 * D)).astype(np.int64)
        _u, first, inverse = np.unique(key, axis=0, return_index=True, return_inverse=True)
        inverse = inverse.reshape(-1)
        stations = centre[sel][first]
        Tb_station = np.full(len(stations), np.nan)
        for s_i, p in enumerate(stations):
            idx = np.array(tree.query_ball_point(p, 1.6 * D))
            if idx.size < 4:
                continue
            nn = nodes[idx]
            w = node_vol[nn]
            e = (w[:, None] * U[nn]).sum(axis=0)
            speed = np.linalg.norm(e)
            if speed <= 0.0:
                continue
            e /= speed
            slab = np.abs((X[nn] - p) @ e) < 0.15 * D
            un = U[nn[slab]] @ e
            den = (w[slab] * un).sum()
            if den > 0.0:
                Tb_station[s_i] = (w[slab] * un * T[nn[slab]]).sum() / den
        T_b[sel] = Tb_station[inverse]

    with np.errstate(divide="ignore", invalid="ignore"):
        h_in = q_s / (T_w - args.t_in)
        h_b = q_s / (T_w - T_b)
    Nu_in = h_in * dh / args.k_fluid
    Nu_b = h_b * dh / args.k_fluid

    # ================================
    # CHECKS AND PER-CHANNEL SUMMARY
    # ================================
    Q_s = (q_s * area).sum()
    Q_f = (q_f * area).sum()
    qin = args.qin
    if qin is None and args.log:
        qin = heat_input_from_log(args.log)
    print("HEAT INTO THE COOLANT: consistent %.6g W | element gradients: steel %.6g W, coolant %.6g W"
          % (Q_s, (q_sg * area).sum(), Q_f))
    print("STEEL NODES OFF THE INTERFACE: sum |nodal residual| = %.3g W (0 AT STEADY STATE)" % resid)
    if qin:
        print("HEAT INPUT           : %.6g W x f %.4g -> consistent / input = %.4f (1 AT STEADY STATE)"
              % (qin, args.heat_scale, Q_s / (qin * args.heat_scale)))
    rows = []
    print("%-8s %10s %10s %10s %10s %10s %10s" % ("channel", "area m2", "Q W", "h_in", "h_b", "Nu_in", "Nu_b"))
    for k in range(1, nch + 1):
        m = ch == k
        good = m & np.isfinite(h_b) & (np.abs(T_w - T_b) > 0.05)
        A = area[m].sum()
        Qk = (q_s[m] * area[m]).sum()
        # AREA-WEIGHTED MEANS; h FROM THE MEAN FLUX OVER THE MEAN DIFFERENCE
        hin = Qk / ((T_w[m] - args.t_in) * area[m]).sum()
        hb = (q_s[good] * area[good]).sum() / ((T_w[good] - T_b[good]) * area[good]).sum() if good.any() else np.nan
        row = dict(channel=k, D_h=args.dh[k - 1], area_m2=A, Q_W=Qk, h_in=hin, h_b=hb,
                   Nu_in=hin * args.dh[k - 1] / args.k_fluid, Nu_b=hb * args.dh[k - 1] / args.k_fluid)
        rows.append(row)
        print("%-8d %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g"
              % (k, A, Qk, row["h_in"], row["h_b"], row["Nu_in"], row["Nu_b"]))
    with open(out + "_channels.csv", "w") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ================================
    # SURFACE FILE FOR PARAVIEW
    # ================================
    nodes, inv = np.unique(tri.reshape(-1), return_inverse=True)
    cells = inv.reshape(-1, 3)
    fields = dict(q_w=q_s, q_w_solid_gradient=q_sg, q_w_fluid_gradient=q_f, T_w=T_w, T_b=T_b, h_in=h_in, h_b=h_b, Nu_in=Nu_in,
                  Nu_b=Nu_b, channel=ch.astype(float))
    with open(out + "_interface.vtu", "w") as v:
        v.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n'
                '<UnstructuredGrid><Piece NumberOfPoints="%d" NumberOfCells="%d">\n' % (len(nodes), nf))
        v.write('<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        np.savetxt(v, X[nodes], fmt="%.9g")
        v.write("</DataArray></Points>\n<Cells>\n")
        v.write('<DataArray type="Int64" Name="connectivity" format="ascii">\n')
        np.savetxt(v, cells, fmt="%d")
        v.write('</DataArray>\n<DataArray type="Int64" Name="offsets" format="ascii">\n')
        np.savetxt(v, 3 * np.arange(1, nf + 1), fmt="%d")
        v.write('</DataArray>\n<DataArray type="UInt8" Name="types" format="ascii">\n')
        np.savetxt(v, np.full(nf, 5), fmt="%d")
        v.write('</DataArray>\n</Cells>\n<CellData Scalars="h_b">\n')
        for name, val in fields.items():
            v.write('<DataArray type="Float64" Name="%s" format="ascii">\n' % name)
            np.savetxt(v, np.nan_to_num(val, nan=0.0, posinf=0.0, neginf=0.0), fmt="%.7g")
            v.write("</DataArray>\n")
        v.write("</CellData>\n</Piece></UnstructuredGrid></VTKFile>\n")
    print("WROTE  : %s_interface.vtu, %s_channels.csv" % (out, out))

    # ================================
    # PROFILES ALONG THE HEATED WALL OF EACH CHANNEL: THE INTERFACE FACES
    # WHOSE NORMAL (SOLID -> FLUID) POINTS AWAY FROM THE FIRST WALL (-y),
    # I.E. THE CHANNEL CEILING FACING THE HEAT SOURCE
    # ================================
    ax_i = "xyz".index(args.axis)
    prof = []
    for k in range(1, nch + 1):
        m = (ch == k) & (n[:, 1] < -0.9)
        if not m.any():
            continue
        s = centre[m, ax_i]
        D = args.dh[k - 1]
        width = args.bin * D
        edges = np.arange(s.min(), s.max() + width, width)
        bins = np.digitize(s, edges)
        for bidx in np.unique(bins):
            mb = bins == bidx
            Ab = area[m][mb]
            ok = np.isfinite(T_b[m][mb])
            if not ok.any():
                continue
            # BIN VALUES: HEAT / AREA AND h = SUM(q A) / SUM((T_w - T_b) A)
            Q_bin = (q_s[m][mb] * Ab).sum()
            Tw_bin = (T_w[m][mb] * Ab).sum() / Ab.sum()
            Tb_bin = (T_b[m][mb][ok] * Ab[ok]).sum() / Ab[ok].sum()
            hb_bin = (q_s[m][mb][ok] * Ab[ok]).sum() / ((T_w[m][mb][ok] - T_b[m][mb][ok]) * Ab[ok]).sum()
            prof.append(dict(channel=k, s=float((s[mb] * Ab).sum() / Ab.sum()), q_w=float(Q_bin / Ab.sum()),
                             T_w=float(Tw_bin), T_b=float(Tb_bin), h_b=float(hb_bin),
                             Nu_b=float(hb_bin * D / args.k_fluid)))
    if prof:
        with open(out + "_profiles.csv", "w") as handle:
            w = csv.DictWriter(handle, fieldnames=list(prof[0].keys()))
            w.writeheader()
            w.writerows(prof)
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
            fig, axs = plt.subplots(2, 2, figsize=(10, 6.5), sharex=True)
            colours = ["#1b6ca8", "#d1495b", "#2a9d8f", "#e9a03b"]
            for k in range(1, nch + 1):
                P = [p for p in prof if p["channel"] == k]
                if not P:
                    continue
                s = np.array([p["s"] for p in P]) * 1e3
                col = colours[(k - 1) % 4]
                axs[0, 0].plot(s, [p["q_w"] / 1e3 for p in P], color=col, label="channel %d" % k)
                axs[0, 1].plot(s, [p["T_w"] for p in P], color=col)
                axs[0, 1].plot(s, [p["T_b"] for p in P], color=col, ls="--")
                axs[1, 0].plot(s, [p["h_b"] for p in P], color=col)
                axs[1, 1].plot(s, [p["Nu_b"] for p in P], color=col)
            axs[0, 0].set_ylabel(r"$q_w$  [kW m$^{-2}$]")
            axs[0, 1].set_ylabel(r"$T$  [K]  (solid: wall $T_w$, dashed: bulk $T_b$)")
            axs[1, 0].set_ylabel(r"$h_b$  [W m$^{-2}$ K$^{-1}$]")
            axs[1, 1].set_ylabel(r"$Nu_b = h_b D_h / k_f$")
            for axx in axs[1]:
                axx.set_xlabel("%s  [mm]" % args.axis)
            axs[0, 0].legend(frameon=False)
            fig.suptitle("Heated wall (channel ceiling facing the first wall) of each coolant channel",
                         fontsize=11)
            fig.tight_layout()
            fig.savefig(out + "_profiles.png", dpi=200)
            print("WROTE  : %s_profiles.png, %s_profiles.csv" % (out, out))
        except ImportError:
            print("NOTE   : matplotlib not available; profiles written as CSV only")


if __name__ == "__main__":
    main()
