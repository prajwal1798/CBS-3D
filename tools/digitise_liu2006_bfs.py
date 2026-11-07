#!/usr/bin/env python3
# ================================
# FILE      : digitise_liu2006_bfs.py
# PURPOSE   : EXACT EXTRACTION OF THE BACKWARD-FACING-STEP REFERENCE DATA FROM
#             THE VECTOR GRAPHICS OF
#
#               P. NITHIARASU, C.-B. LIU, COMPUT. METHODS APPL. MECH. ENGRG.
#               195 (2006) 2961-2982, FIG. 4(b) (SA MODEL, Re = 3025)
#
#             THE FIGURE IS STORED AS PDF PATHS, SO NO PIXEL READING IS
#             INVOLVED: THE EXPERIMENTAL POINTS (DENHAM ET AL. 1975) ARE
#             1.63 pt SQUARES (CENTRE = DATA POINT) AND THE SA SOLUTION IS A
#             POLYLINE. PAGE COORDINATES ARE MAPPED TO DATA COORDINATES WITH
#             THE AXIS FRAME (u_plot = 0..14, y = 0..3).
#
#             EACH PROFILE IS DRAWN SHIFTED: u_plot = u/U0 + OFFSET, OFFSETS
#             0, 2, 4, 6, 8, 12 (DISPLAY SHIFTS, NOT POSITIONS). STATIONS:
#             DENHAM, BRIARD & PATRICK, J. PHYS. E 8 (1975) 681, FIG. 4 HAS
#             x = -20, 0, 15, 30, 60, 90, 120 mm (h = 15 mm). LIU'S SECOND TO
#             SIXTH PROFILES ARE x = 15, 30, 60, 90, 120 mm, I.E. x/h = 1, 2,
#             4, 6, 8: THE DEEPEST REVERSAL IS IN HIS THIRD PROFILE (-0.235)
#             AND DENHAM'S IS AT 30 mm; 15 mm IS WEAKER (-0.12); 60 AND 90 mm
#             ARE ~0 AT THE WALL. READING THE OFFSETS AS x/h WOULD PUT THE
#             DEEPEST REVERSAL AT 60 mm, WHERE DENHAM HAS NONE. THE FIRST
#             PROFILE SPANS y = 1..3 ONLY: x/h = 0 (DENHAM'S -20 mm PROFILE
#             IS NEARLY THE SAME).
#
# OUTPUT    : <out>/liu2006_fig4b_exp.csv   station,x_over_h,y,u   (DENHAM)
#             <out>/liu2006_fig4b_sa.csv    station,x_over_h,y,u   (LIU SA)
#             <out>/inlet_denham_station0.txt  "y u" TABLE FOR BC 507
#             MASS-FLUX CHECK: THE SA FLUX INT u dy IS 1.932-1.955 AT ALL
#             SIX STATIONS (+-0.6 %), AS A CONVERGED SOLUTION REQUIRES.
#
# USAGE     : python3 tools/digitise_liu2006_bfs.py PAPER.pdf --out DIR
#             (NEEDS PyMuPDF: pip install pymupdf)
# ================================

import argparse
import os
import sys

import numpy as np

PAGE = 11                                   # ZERO-BASED: JOURNAL PAGE 2972
OFFSETS = [0.0, 2.0, 4.0, 6.0, 8.0, 12.0]
X_OVER_H = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0]
BREAKS = [1.6, 3.45, 5.5, 7.45, 10.0]       # u_plot BOUNDARIES BETWEEN PROFILES


def frame_b(drawings):
    # PANEL (b): THE AXIS RECTANGLE WHOSE TOP EDGE IS NEAR y = 353.6 pt
    for d in drawings:
        for it in d["items"]:
            if it[0] == "re":
                r = it[1]
                if r.width > 200 and r.height > 100 and 340 < r.y0 < 360:
                    return r
    raise RuntimeError("PANEL (b) FRAME NOT FOUND")


def legend_rows(page):
    # (LABEL, x_right, y_centre) OF THE 'Exp.' AND 'SA model' LEGEND LABELS;
    # EACH SAMPLE (A SQUARE, OR A HORIZONTAL LINE) SITS JUST RIGHT OF IT
    rows = []
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            t = "".join(s["text"] for s in l["spans"]).strip()
            if t in ("Exp.", "SA model") and 340 < l["bbox"][1] < 400:
                x0, y0, x1, y1 = l["bbox"]
                rows.append((t, x1, 0.5 * (y0 + y1)))
    return rows


def is_legend_square(rows, cx, cy):
    return any(t == "Exp." and x1 + 2.0 < cx < x1 + 12.0 and abs(cy - yc) < 1.5 for t, x1, yc in rows)


def is_legend_line(rows, p, q):
    return any(t == "SA model" and abs(p.y - q.y) < 1e-3 and abs(p.y - yc) < 2.0
               and min(p.x, q.x) > x1 and abs(p.x - q.x) > 5.0 for t, x1, yc in rows)


def is_tick(fr, p, q):
    # SHORT AXIS-ALIGNED SEGMENT TOUCHING THE FRAME (TICK MARK)
    if abs(p.x - q.x) > 1e-3 and abs(p.y - q.y) > 1e-3:
        return False
    if abs(p.x - q.x) + abs(p.y - q.y) > 2.5:     # TICKS ARE ~1.6 pt
        return False
    near = lambda v, e: abs(v - e) < 0.3
    return any(near(c.x, fr.x0) or near(c.x, fr.x1) or near(c.y, fr.y0) or near(c.y, fr.y1) for c in (p, q))


def station_of(u_plot):
    return int(np.searchsorted(BREAKS, u_plot))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()

    import pymupdf

    page = pymupdf.open(args.pdf)[PAGE]
    drawings = page.get_drawings()
    fr = frame_b(drawings)
    legend = legend_rows(page)

    def to_data(x, y):
        return 14.0 * (x - fr.x0) / fr.width, 3.0 * (fr.y1 - y) / fr.height

    squares = []
    segments = []
    for d in drawings:
        r = d["rect"]
        if not (fr.x0 - 2 <= r.x0 and r.x1 <= fr.x1 + 2 and fr.y0 - 2 <= r.y0 and r.y1 <= fr.y1 + 2):
            continue
        if r.width > 200 and r.height > 100 and len(d["items"]) <= 2:
            continue                                # FRAME AND GRID LINES
        for it in d["items"]:
            if it[0] == "re" and it[1].width < 3.0:
                cx = 0.5 * (it[1].x0 + it[1].x1)
                cy = 0.5 * (it[1].y0 + it[1].y1)
                if not is_legend_square(legend, cx, cy):
                    squares.append(to_data(cx, cy))
            elif it[0] == "l":
                p, q = it[1], it[2]
                if abs(p.x - q.x) + abs(p.y - q.y) < 0.5:
                    continue                        # MARKER DOTS
                if is_legend_line(legend, p, q) or is_tick(fr, p, q):
                    continue
                if abs(p.x - q.x) < 1e-3 and (abs(p.x - fr.x0) < 0.2 or abs(p.x - fr.x1) < 0.2) \
                        and abs(p.y - q.y) > 6.0:
                    continue                        # FRAME EDGE
                segments.append((to_data(p.x, p.y), to_data(q.x, q.y)))

    # ================================
    # CHAIN THE SA SEGMENTS INTO POLYLINES
    # ================================
    lines = []
    for a, b in segments:
        if lines and np.hypot(lines[-1][-1][0] - a[0], lines[-1][-1][1] - a[1]) < 1e-3:
            lines[-1].append(b)
        else:
            lines.append([a, b])

    sa = []
    for ln in lines:
        arr = np.array(ln)
        k = station_of(np.median(arr[:, 0]))
        for u_plot, y in arr:
            # THE CURVE OF THE FIRST STATION ALSO RUNS DOWN THE u = 0 LINE
            # BELOW THE STEP (y < 1); KEEP ONLY THE FLOW PART
            if k == 0 and y < 1.0 - 1e-6:
                continue
            sa.append((k, y, u_plot - OFFSETS[k]))

    exp = [(station_of(u), y, u - OFFSETS[station_of(u)]) for u, y in squares]

    os.makedirs(args.out, exist_ok=True)
    src = "Nithiarasu & Liu, CMAME 195 (2006) 2961, Fig. 4(b); vector extraction by tools/digitise_liu2006_bfs.py"
    for name, rows, what in (("liu2006_fig4b_exp.csv", exp, "Denham, Briard & Patrick (1975) data as plotted"),
                             ("liu2006_fig4b_sa.csv", sa, "Liu SA solution")):
        rows = sorted(rows, key=lambda r: (r[0], r[1]))
        with open(os.path.join(args.out, name), "w") as f:
            f.write("# %s\n# %s\n# u = u/U0 (U0 = MEAN UPSTREAM VELOCITY), y/h FROM THE DOWNSTREAM FLOOR\n" % (what, src))
            f.write("station,x_over_h,y,u\n")
            for k, y, u in rows:
                f.write("%d,%.1f,%.5f,%.5f\n" % (k, X_OVER_H[k], y, u))
        counts = [sum(1 for r in rows if r[0] == k) for k in range(len(OFFSETS))]
        print("%-24s %4d points, per station %s" % (name, len(rows), counts), file=sys.stderr)

    # ================================
    # INLET TABLE FOR BC 507 (CBS3D_BFS_INLET_PROFILE): THE MEASURED PROFILE
    # OF THE FIRST STATION, y = 1..3, WITH u = 0 AT BOTH WALLS. POINTS BELOW
    # THE UPSTREAM FLOOR (PLOTTING ARTEFACTS) AND THE PLOTTED WALL POINTS ARE
    # REPLACED BY THE EXACT WALL VALUES.
    # ================================
    first = sorted((y, u) for k, y, u in exp if k == 0 and 1.0 + 1e-3 < y < 3.0 - 1e-3)
    table = [(1.0, 0.0)] + first + [(3.0, 0.0)]
    ys = np.array([r[0] for r in table])
    us = np.array([r[1] for r in table])
    with open(os.path.join(args.out, "inlet_denham_station0.txt"), "w") as f:
        f.write("# BC 507 INLET TABLE: y/h  u/U0  (LINEAR INTERPOLATION; 0 OUTSIDE)\n")
        f.write("# Denham, Briard & Patrick (1975), first station, as plotted in %s\n" % src)
        f.write("# FLUX INT u dy = %.5f (MEAN u = %.5f OVER y = 1..3)\n"
                % (np.trapezoid(us, ys), np.trapezoid(us, ys) / 2.0))
        for y, u in table:
            f.write("%.5f %.5f\n" % (y, u))
    print("inlet_denham_station0.txt %d points, flux %.4f" % (len(table), np.trapezoid(us, ys)), file=sys.stderr)


if __name__ == "__main__":
    main()
