#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CBS3D SA flat-plate publication post-processing.

Python 3.6 compatible.

This script reads the CSV output already produced by
postprocess_flatplate_sa_pvtu_py36.py, downloads the official NASA/TMBWG
flat-plate SA reference files, and generates publication-style comparison
figures.

Required CBS3D files in --post-dir:
  flatplate_cf_binned.csv
  profile_x*.csv

Optional, auto-detected one directory above --post-dir:
  flatplate_distributed_residuals.csv

Official NASA/TMBWG files downloaded automatically:
  FlatPlate/SA/cf_plate.dat
  FlatPlate/SA/cf_plate_fun3d_tri.dat
  FlatPlate/SA/cf_convergence.dat
  FlatPlate/SA/cf_convergence_fun3d_tri.dat
  FlatPlate_validation/cf_incomp_results_sa.dat
  FlatPlate/SA/flatplate_u+y+.dat
  FlatPlate/SA/u+y+theory.dat
  FlatPlate/SA/flatplate_u.dat
  FlatPlate/SA/mut_0.97.dat

Outputs:
  paper_plots/*.png
  paper_plots/*.pdf
  paper_plots/validation_summary.txt
  paper_plots/validation_summary.csv
"""

from __future__ import print_function

import argparse
import csv
import glob
import math
import os
import re
import shutil
import subprocess
import sys

try:
    from urllib.request import Request, urlopen
except ImportError:
    from urllib2 import Request, urlopen

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TMR_BASE_RAW = "https://raw.githubusercontent.com/TMBWG/turbmodels/main/"
TMR_BASE_WEB = "https://tmbwg.github.io/turbmodels/"

TMR_FILES = {
    "cf_plate.dat": "FlatPlate/SA/cf_plate.dat",
    "cf_plate_fun3d_tri.dat": "FlatPlate/SA/cf_plate_fun3d_tri.dat",
    "cf_convergence.dat": "FlatPlate/SA/cf_convergence.dat",
    "cf_convergence_fun3d_tri.dat": "FlatPlate/SA/cf_convergence_fun3d_tri.dat",
    "cf_incomp_results_sa.dat": "FlatPlate_validation/cf_incomp_results_sa.dat",
    "flatplate_u+y+.dat": "FlatPlate/SA/flatplate_u+y+.dat",
    "u+y+theory.dat": "FlatPlate/SA/u+y+theory.dat",
    "flatplate_u.dat": "FlatPlate/SA/flatplate_u.dat",
    "mut_0.97.dat": "FlatPlate/SA/mut_0.97.dat",
}

TARGET_X1 = 0.97008
TARGET_X2 = 1.90334
CBS_EQUIV_N = 13056.0  # (137-1)*(97-1)
CBS_EQUIV_H = math.sqrt(1.0 / CBS_EQUIV_N)


# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------

def apply_paper_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9.5,
        "mathtext.fontset": "cm",
        "axes.labelsize": 10.5,
        "legend.fontsize": 8.5,
        "xtick.labelsize": 9.0,
        "ytick.labelsize": 9.0,
        "axes.linewidth": 0.75,
        "lines.linewidth": 1.25,
        "lines.markersize": 4.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "savefig.dpi": 400,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def style_axes(ax, grid=False):
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    if grid:
        ax.grid(True, which="major", linestyle=":", linewidth=0.45, alpha=0.45)


def save_both(fig, out_dir, stem):
    png = os.path.join(out_dir, stem + ".png")
    pdf = os.path.join(out_dir, stem + ".pdf")
    fig.savefig(png, dpi=400, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


# ---------------------------------------------------------------------------
# File utilities
# ---------------------------------------------------------------------------

def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def download_one(url, dst):
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    response = urlopen(req, timeout=40)
    data = response.read()
    if len(data) < 20:
        raise RuntimeError("downloaded file is unexpectedly small")
    with open(dst, "wb") as f:
        f.write(data)


def curl_download(url, dst):
    curl = shutil.which("curl")
    if curl is None:
        raise RuntimeError("curl is not available")
    cmd = [curl, "-L", "--fail", "--silent", "--show-error", url, "-o", dst]
    subprocess.check_call(cmd)
    if (not os.path.isfile(dst)) or os.path.getsize(dst) < 20:
        raise RuntimeError("curl output is missing or too small")


def ensure_tmr_references(ref_dir):
    ensure_dir(ref_dir)
    paths = {}

    for local_name in sorted(TMR_FILES.keys()):
        rel = TMR_FILES[local_name]
        dst = os.path.join(ref_dir, local_name)
        paths[local_name] = dst

        if os.path.isfile(dst) and os.path.getsize(dst) > 20:
            print("[ref] present  : {}".format(local_name))
            continue

        urls = [
            TMR_BASE_RAW + rel,
            TMR_BASE_WEB + rel,
        ]

        ok = False
        last_error = None

        for url in urls:
            try:
                print("[ref] download : {}".format(url))
                download_one(url, dst)
                ok = True
                break
            except Exception as exc:
                last_error = exc
                try:
                    curl_download(url, dst)
                    ok = True
                    break
                except Exception as exc2:
                    last_error = exc2

        if not ok:
            raise RuntimeError(
                "Could not download official TMR file '{}': {}".format(
                    local_name, last_error
                )
            )

    return paths


# ---------------------------------------------------------------------------
# Tecplot ASCII parser
# ---------------------------------------------------------------------------

def zone_title(line):
    m = re.search(r't\s*=\s*"([^"]+)"', line, flags=re.IGNORECASE)
    if m:
        return m.group(1).strip()
    return "zone"


def parse_tecplot(path):
    variables = []
    zones = []
    current_title = "default"
    current_rows = []

    def flush():
        if current_rows:
            arr = np.asarray(current_rows, dtype=float)
            zones.append({
                "title": current_title,
                "data": arr,
            })

    with open(path, "r") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#"):
                continue

            low = s.lower()

            if low.startswith("variables"):
                variables = re.findall(r'"([^"]+)"', s)
                continue

            if low.startswith("zone"):
                flush()
                current_title = zone_title(s)
                current_rows = []
                continue

            parts = s.replace(",", " ").split()
            vals = []
            numeric = True
            for token in parts:
                try:
                    vals.append(float(token))
                except ValueError:
                    numeric = False
                    break

            if numeric and vals:
                if variables and len(vals) < len(variables):
                    continue
                if variables:
                    vals = vals[:len(variables)]
                current_rows.append(vals)

    flush()

    if not zones:
        raise RuntimeError("No numeric Tecplot zones parsed from {}".format(path))

    return {
        "variables": variables,
        "zones": zones,
    }


def find_zone(zones, contains):
    key = contains.lower()
    for z in zones:
        if key in z["title"].lower():
            return z
    return None


def find_station_zone(zones, target_x):
    best = None
    best_dist = None
    for z in zones:
        nums = re.findall(r'[-+]?(?:\d+\.\d+|\d+)', z["title"])
        if not nums:
            continue
        try:
            x = float(nums[-1])
        except ValueError:
            continue
        d = abs(x - target_x)
        if best is None or d < best_dist:
            best = z
            best_dist = d
    return best


# ---------------------------------------------------------------------------
# CBS3D CSV readers
# ---------------------------------------------------------------------------

def read_csv_rows(path):
    with open(path, "r") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError("CSV has no data: {}".format(path))
    return rows


def column(rows, name):
    return np.asarray([float(r[name]) for r in rows], dtype=float)


def profile_x_from_filename(path):
    name = os.path.basename(path)
    m = re.search(r"profile_x([0-9mp]+)\.csv$", name)
    if not m:
        return math.nan
    return float(m.group(1).replace("m", "-").replace("p", "."))


def load_cbs_profiles(post_dir):
    out = []
    files = sorted(glob.glob(os.path.join(post_dir, "profile_x*.csv")))
    if not files:
        raise RuntimeError("No profile_x*.csv found in {}".format(post_dir))

    for path in files:
        rows = read_csv_rows(path)
        rec = {
            "path": path,
            "x": profile_x_from_filename(path),
            "y": column(rows, "y"),
            "u": column(rows, "u"),
            "v": column(rows, "v"),
            "w": column(rows, "w"),
            "mut": column(rows, "mu_t_over_mu"),
            "wall_distance": column(rows, "wall_distance"),
            "y_plus": column(rows, "y_plus"),
            "u_plus": column(rows, "u_plus"),
        }
        out.append(rec)

    return out


def nearest_cbs_profile(profiles, target_x):
    return min(profiles, key=lambda p: abs(p["x"] - target_x))


def load_cbs_cf(post_dir):
    path = os.path.join(post_dir, "flatplate_cf_binned.csv")
    rows = read_csv_rows(path)
    x = column(rows, "x")
    cf = column(rows, "Cf")
    order = np.argsort(x)
    return x[order], cf[order]


def read_summary(post_dir):
    path = os.path.join(post_dir, "summary.txt")
    d = {}
    if not os.path.isfile(path):
        return d

    with open(path, "r") as f:
        for raw in f:
            s = raw.strip()
            if " = " not in s:
                continue
            key, val = s.split(" = ", 1)
            d[key.strip()] = val.strip()
    return d


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------

def real_profile_mask(profile):
    return (
        np.isfinite(profile["y"]) &
        np.isfinite(profile["y_plus"]) &
        np.isfinite(profile["u_plus"]) &
        (profile["y"] > 1.0e-12) &
        (profile["y_plus"] > 0.0)
    )


def first_real_yplus(profile):
    mask = real_profile_mask(profile)
    if not np.any(mask):
        return math.nan
    y = profile["y"][mask]
    yp = profile["y_plus"][mask]
    idx = np.argmin(y)
    return float(yp[idx])


def delta99(y, u_over_uinf):
    mask = np.isfinite(y) & np.isfinite(u_over_uinf) & (y >= 0.0)
    yy = y[mask]
    uu = u_over_uinf[mask]
    if yy.size < 3:
        return math.nan

    order = np.argsort(yy)
    yy = yy[order]
    uu = uu[order]

    n_outer = max(3, int(math.ceil(0.10 * yy.size)))
    ue = float(np.median(uu[-n_outer:]))
    if (not math.isfinite(ue)) or ue <= 0.0:
        return math.nan

    target = 0.99 * ue
    hit = np.where(uu >= target)[0]
    if hit.size == 0:
        return math.nan

    i = int(hit[0])
    if i == 0:
        return float(yy[0])

    y0 = float(yy[i - 1])
    y1 = float(yy[i])
    u0 = float(uu[i - 1])
    u1 = float(uu[i])

    if abs(u1 - u0) < 1.0e-15:
        return y1

    a = (target - u0) / (u1 - u0)
    return y0 + a * (y1 - y0)


def interp_at(x, y, xq):
    mask = np.isfinite(x) & np.isfinite(y)
    xx = x[mask]
    yy = y[mask]
    order = np.argsort(xx)
    return float(np.interp(xq, xx[order], yy[order]))


# ---------------------------------------------------------------------------
# Reference extraction
# ---------------------------------------------------------------------------

def ref_cf_distribution(paths):
    q = parse_tecplot(paths["cf_plate.dat"])
    tri = parse_tecplot(paths["cf_plate_fun3d_tri.dat"])

    result = []

    z = find_zone(q["zones"], "CFL3D")
    if z is not None:
        result.append(("TMR CFL3D (M=0.2)", z["data"][:, 0], z["data"][:, 1]))

    z = find_zone(q["zones"], "FUN3D")
    if z is not None:
        result.append(("TMR FUN3D quad (M=0.2)", z["data"][:, 0], z["data"][:, 1]))

    for z in tri["zones"]:
        result.append(("TMR FUN3D triangles (M=0.2)", z["data"][:, 0], z["data"][:, 1]))
        break

    return result


def ref_cf_convergence(paths):
    result = []

    q = parse_tecplot(paths["cf_convergence.dat"])
    for name in ["CFL3D", "FUN3D"]:
        z = find_zone(q["zones"], name)
        if z is not None and z["data"].shape[1] >= 4:
            result.append((
                "TMR {} compressible".format(name),
                z["data"][:, 2],
                z["data"][:, 3],
            ))

    tri = parse_tecplot(paths["cf_convergence_fun3d_tri.dat"])
    if tri["zones"] and tri["zones"][0]["data"].shape[1] >= 4:
        z = tri["zones"][0]
        result.append((
            "TMR FUN3D triangles compressible",
            z["data"][:, 2],
            z["data"][:, 3],
        ))

    inc = parse_tecplot(paths["cf_incomp_results_sa.dat"])
    for z in inc["zones"]:
        if z["data"].shape[1] >= 4:
            result.append((
                "TMR " + z["title"],
                z["data"][:, 2],
                z["data"][:, 3],
            ))

    return result


def ref_incompressible_same_grid(paths):
    inc = parse_tecplot(paths["cf_incomp_results_sa.dat"])
    out = {}
    for z in inc["zones"]:
        data = z["data"]
        if data.shape[1] < 4:
            continue
        idx = int(np.argmin(np.abs(data[:, 0] - CBS_EQUIV_N)))
        out[z["title"]] = {
            "N": float(data[idx, 0]),
            "h": float(data[idx, 2]),
            "cf": float(data[idx, 3]),
        }
    return out


def ref_uplus(paths, target_x):
    q = parse_tecplot(paths["flatplate_u+y+.dat"])
    z = find_station_zone(q["zones"], target_x)
    if z is None:
        raise RuntimeError("No TMR u+/y+ station found near x={}".format(target_x))
    logyp = z["data"][:, 0]
    up = z["data"][:, 1]
    return z["title"], np.power(10.0, logyp), up


def ref_theory(paths):
    q = parse_tecplot(paths["u+y+theory.dat"])
    z = q["zones"][0]
    d = z["data"]
    if d.shape[1] < 5:
        raise RuntimeError("Unexpected u+y+theory.dat format")
    up = d[:, 0]
    return {
        "u_plus": up,
        "yplus_inner": np.power(10.0, d[:, 1]),
        "yplus_log": np.power(10.0, d[:, 2]),
        "yplus_spalding": np.power(10.0, d[:, 4]),
    }


def ref_velocity(paths, target_x):
    q = parse_tecplot(paths["flatplate_u.dat"])
    z = find_station_zone(q["zones"], target_x)
    if z is None:
        raise RuntimeError("No TMR velocity station found near x={}".format(target_x))
    return z["title"], z["data"][:, 0], z["data"][:, 1]


def ref_mut(paths):
    q = parse_tecplot(paths["mut_0.97.dat"])
    out = []
    for z in q["zones"]:
        d = z["data"]
        if d.shape[1] >= 3:
            out.append((z["title"], d[:, 1], d[:, 2]))
    return out


# ---------------------------------------------------------------------------
# Plot functions
# ---------------------------------------------------------------------------

def plot_cf_distribution(cbs_x, cbs_cf, refs, out_dir):
    styles = ["--", "-.", ":"]
    markers = [None, None, None]

    # Main engineering-range view
    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    ax.plot(cbs_x, cbs_cf, "-", label="CBS3D", linewidth=1.45)

    for i, item in enumerate(refs):
        ax.plot(
            item[1], item[2],
            linestyle=styles[i % len(styles)],
            marker=markers[i % len(markers)],
            linewidth=1.05,
            label=item[0],
        )

    ax.axvline(TARGET_X1, linestyle=":", linewidth=0.7)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$C_f$")
    ax.set_xlim(0.05, 1.95)

    mask = (cbs_x >= 0.05) & (cbs_x <= 1.95)
    vals = list(cbs_cf[mask])
    for item in refs:
        m = (item[1] >= 0.05) & (item[1] <= 1.95)
        vals.extend(list(item[2][m]))
    vals = np.asarray([v for v in vals if math.isfinite(float(v))], dtype=float)
    if vals.size:
        lo = float(np.min(vals))
        hi = float(np.max(vals))
        pad = 0.08 * max(hi - lo, 1.0e-6)
        ax.set_ylim(lo - pad, hi + pad)

    style_axes(ax, grid=False)
    ax.legend(loc="best")
    fig.tight_layout()
    save_both(fig, out_dir, "01_cf_distribution_main")

    # Full plate diagnostic
    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    ax.plot(cbs_x, cbs_cf, "-", label="CBS3D", linewidth=1.35)
    for i, item in enumerate(refs):
        ax.plot(
            item[1], item[2],
            linestyle=styles[i % len(styles)],
            linewidth=1.0,
            label=item[0],
        )
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$C_f$")
    ax.set_xlim(0.0, 2.0)
    style_axes(ax, grid=False)
    ax.legend(loc="best")
    fig.tight_layout()
    save_both(fig, out_dir, "01b_cf_distribution_full")


def plot_cf_convergence(cbs_cf_x097, refs, out_dir):
    styles = ["-", "--", "-.", ":", "--"]
    markers = ["o", "s", "^", "D", "v"]

    fig, ax = plt.subplots(figsize=(5.6, 3.9))

    for i, item in enumerate(refs):
        h = np.asarray(item[1], dtype=float)
        cf = np.asarray(item[2], dtype=float)
        order = np.argsort(h)
        ax.plot(
            h[order], cf[order],
            linestyle=styles[i % len(styles)],
            marker=markers[i % len(markers)],
            linewidth=1.0,
            markersize=3.6,
            label=item[0],
        )

    ax.plot(
        [CBS_EQUIV_H], [cbs_cf_x097],
        marker="*", linestyle="None", markersize=7.0,
        label="CBS3D coarse (137x97-derived TET4)",
    )

    ax.set_xlabel(r"$h=\sqrt{1/N}$")
    ax.set_ylabel(r"$C_f(x=0.97008)$")
    ax.set_xlim(left=0.0)
    style_axes(ax, grid=False)
    ax.legend(loc="best")
    fig.tight_layout()
    save_both(fig, out_dir, "02_cf_x097_grid_convergence")


def plot_uplus(profile, target_x, paths, out_dir):
    zone_name, tmr_yp, tmr_up = ref_uplus(paths, target_x)
    theory = ref_theory(paths)

    mask = real_profile_mask(profile)
    cbs_yp = profile["y_plus"][mask]
    cbs_up = profile["u_plus"][mask]
    order = np.argsort(cbs_yp)
    cbs_yp = cbs_yp[order]
    cbs_up = cbs_up[order]

    fig, ax = plt.subplots(figsize=(4.8, 3.9))

    ax.semilogx(
        cbs_yp, cbs_up, "-",
        linewidth=1.45,
        label="CBS3D (x={:.5f})".format(profile["x"]),
    )
    ax.semilogx(
        tmr_yp, tmr_up, "--",
        linewidth=1.1,
        label="TMR CFL3D ({})".format(re.sub(r"\s+", " ", zone_name)),
    )

    inner_mask = (
        np.isfinite(theory["yplus_inner"]) &
        np.isfinite(theory["u_plus"]) &
        (theory["yplus_inner"] > 0.0)
    )
    log_mask = (
        np.isfinite(theory["yplus_log"]) &
        np.isfinite(theory["u_plus"]) &
        (theory["yplus_log"] > 0.0)
    )

    ax.semilogx(
        theory["yplus_inner"][inner_mask],
        theory["u_plus"][inner_mask],
        ":",
        linewidth=0.9,
        label=r"TMR inner law",
    )
    ax.semilogx(
        theory["yplus_log"][log_mask],
        theory["u_plus"][log_mask],
        "-.",
        linewidth=0.9,
        label=r"TMR log law ($\kappa=0.41,\ B=5$)",
    )

    ax.set_xlabel(r"$y^+$")
    ax.set_ylabel(r"$u^+$")
    style_axes(ax, grid=False)
    ax.legend(loc="best")
    fig.tight_layout()

    stem = "03_uplus_yplus_x0p97008" if abs(target_x - TARGET_X1) < 0.01 else "04_uplus_yplus_x1p90334"
    save_both(fig, out_dir, stem)


def plot_velocity(profile, target_x, paths, out_dir, uinf):
    zone_name, tmr_u, tmr_y = ref_velocity(paths, target_x)

    cbs_u = profile["u"] / uinf
    cbs_y = profile["y"]

    cbs_d99 = delta99(cbs_y, cbs_u)
    tmr_d99 = delta99(tmr_y, tmr_u)

    candidates = [v for v in [cbs_d99, tmr_d99] if math.isfinite(v) and v > 0.0]
    if candidates:
        y_max = 1.25 * max(candidates)
    else:
        y_max = min(0.12, max(float(np.max(cbs_y)), float(np.max(tmr_y))))

    fig, ax = plt.subplots(figsize=(4.5, 3.9))

    cmask = np.isfinite(cbs_y) & np.isfinite(cbs_u) & (cbs_y >= 0.0) & (cbs_y <= y_max)
    tmask = np.isfinite(tmr_y) & np.isfinite(tmr_u) & (tmr_y >= 0.0) & (tmr_y <= y_max)

    ax.plot(
        cbs_u[cmask], cbs_y[cmask], "-",
        linewidth=1.45,
        label="CBS3D (x={:.5f})".format(profile["x"]),
    )
    ax.plot(
        tmr_u[tmask], tmr_y[tmask], "--",
        linewidth=1.1,
        label="TMR CFL3D ({})".format(re.sub(r"\s+", " ", zone_name)),
    )

    ax.set_xlabel(r"$u/U_\infty$")
    ax.set_ylabel(r"$y$")
    ax.set_xlim(0.0, 1.06)
    ax.set_ylim(0.0, y_max)
    style_axes(ax, grid=False)
    ax.legend(loc="best")
    fig.tight_layout()

    stem = "05_velocity_profile_x0p97008" if abs(target_x - TARGET_X1) < 0.01 else "06_velocity_profile_x1p90334"
    save_both(fig, out_dir, stem)


def plot_mut(profile, paths, out_dir, uinf):
    refs = ref_mut(paths)

    cbs_u = profile["u"] / uinf
    d99 = delta99(profile["y"], cbs_u)
    if math.isfinite(d99) and d99 > 0.0:
        y_max = 1.25 * d99
    else:
        y_max = min(0.12, float(np.max(profile["y"])))

    fig, ax = plt.subplots(figsize=(4.7, 3.9))

    cmask = (
        np.isfinite(profile["mut"]) &
        np.isfinite(profile["y"]) &
        (profile["y"] >= 0.0) &
        (profile["y"] <= y_max)
    )
    ax.plot(
        profile["mut"][cmask], profile["y"][cmask], "-",
        linewidth=1.45,
        label="CBS3D (x={:.5f})".format(profile["x"]),
    )

    styles = ["--", "-."]
    for i, item in enumerate(refs):
        title, y, mut = item
        mask = np.isfinite(y) & np.isfinite(mut) & (y >= 0.0) & (y <= y_max)
        ax.plot(
            mut[mask], y[mask],
            linestyle=styles[i % len(styles)],
            linewidth=1.05,
            label="TMR {}".format(title),
        )

    ax.set_xlabel(r"$\mu_t/\mu_\infty$")
    ax.set_ylabel(r"$y$")
    ax.set_ylim(0.0, y_max)
    ax.set_xlim(left=0.0)
    style_axes(ax, grid=False)
    ax.legend(loc="best")
    fig.tight_layout()
    save_both(fig, out_dir, "07_mut_profile_x0p97")


def find_residual_csv(post_dir):
    candidates = [
        os.path.join(os.path.dirname(post_dir), "flatplate_distributed_residuals.csv"),
        os.path.join(post_dir, "flatplate_distributed_residuals.csv"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    return None


def plot_residual_history(post_dir, out_dir):
    path = find_residual_csv(post_dir)
    if path is None:
        return

    rows = read_csv_rows(path)
    names = rows[0].keys()

    if "iteration" not in names:
        return

    it = column(rows, "iteration")

    series = []
    for key, label in [
        ("u_rel", r"$u$"),
        ("v_rel", r"$v$"),
        ("w_rel", r"$w$"),
        ("p_rel", r"$p$"),
        ("sa_rel", "SA"),
    ]:
        if key in names:
            series.append((key, label))

    if series:
        fig, ax = plt.subplots(figsize=(5.6, 3.8))
        styles = ["-", "--", "-.", ":", "-"]
        for i, item in enumerate(series):
            vals = column(rows, item[0])
            mask = np.isfinite(vals) & (vals > 0.0)
            ax.semilogy(
                it[mask], vals[mask],
                linestyle=styles[i % len(styles)],
                linewidth=1.0,
                label=item[1],
            )
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Relative residual")
        style_axes(ax, grid=False)
        ax.legend(loc="best", ncol=2)
        fig.tight_layout()
        save_both(fig, out_dir, "08_residual_history")


# ---------------------------------------------------------------------------
# Validation summary
# ---------------------------------------------------------------------------

def percent_error(value, ref):
    return 100.0 * (value - ref) / ref


def write_validation_summary(out_dir, post_dir, cbs_cf_x097, profiles, paths):
    same = ref_incompressible_same_grid(paths)
    p1 = nearest_cbs_profile(profiles, TARGET_X1)
    p2 = nearest_cbs_profile(profiles, TARGET_X2)

    rows = []
    rows.append(("CBS3D_Cf_x0p97008", cbs_cf_x097, ""))
    rows.append(("CBS3D_first_real_yplus_x0p97008", first_real_yplus(p1), ""))
    rows.append(("CBS3D_first_real_yplus_x1p90334_nearest", first_real_yplus(p2), ""))
    rows.append(("CBS3D_mesh_x_near_0p97008", p1["x"], ""))
    rows.append(("CBS3D_mesh_x_near_1p90334", p2["x"], ""))

    summary = read_summary(post_dir)
    if "max_abs_w_over_all_piece_nodes" in summary:
        rows.append((
            "CBS3D_max_abs_w",
            float(summary["max_abs_w_over_all_piece_nodes"]),
            "",
        ))
    if "max_mu_t_over_mu_owned_nodes" in summary:
        rows.append((
            "CBS3D_max_mu_t_over_mu",
            float(summary["max_mu_t_over_mu_owned_nodes"]),
            "",
        ))

    for title in sorted(same.keys()):
        ref = same[title]
        err = percent_error(cbs_cf_x097, ref["cf"])
        rows.append((
            "TMR_same_grid_{}_Cf".format(title.replace(" ", "_")),
            ref["cf"],
            "",
        ))
        rows.append((
            "CBS3D_error_vs_{}_percent".format(title.replace(" ", "_")),
            err,
            "%",
        ))

    txt_path = os.path.join(out_dir, "validation_summary.txt")
    with open(txt_path, "w") as f:
        f.write("CBS3D SA flat-plate validation summary\n")
        f.write("=" * 72 + "\n")
        f.write("CBS equivalent TMR logical-grid N = {:.0f}\n".format(CBS_EQUIV_N))
        f.write("CBS equivalent h = {:.12g}\n".format(CBS_EQUIV_H))
        f.write("\n")
        for key, value, unit in rows:
            if unit:
                f.write("{} = {:.12g} {}\n".format(key, value, unit))
            else:
                f.write("{} = {:.12g}\n".format(key, value))
        f.write("\n")
        f.write("Reference provenance: NASA/TMBWG official turbmodels repository.\n")
        f.write("Note: full Cf(x), u+, U-profile, and mu_t/mu profile TMR curves are\n")
        f.write("compressible M=0.2 reference results unless explicitly labelled\n")
        f.write("incompressible. The principal exact incompressible comparison used\n")
        f.write("here is Cf at x=0.97008 from cf_incomp_results_sa.dat.\n")

    csv_path = os.path.join(out_dir, "validation_summary.csv")
    with open(csv_path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["quantity", "value", "unit"])
        for row in rows:
            writer.writerow(row)


# ---------------------------------------------------------------------------
# Internal checks
# ---------------------------------------------------------------------------

def sanity_check_references(paths):
    checks = [
        ("cf_plate.dat", 2),
        ("flatplate_u+y+.dat", 2),
        ("flatplate_u.dat", 2),
        ("mut_0.97.dat", 2),
        ("cf_incomp_results_sa.dat", 2),
    ]

    for name, min_zones in checks:
        parsed = parse_tecplot(paths[name])
        if len(parsed["zones"]) < min_zones:
            raise RuntimeError(
                "{} parsed only {} zone(s), expected at least {}".format(
                    name, len(parsed["zones"]), min_zones
                )
            )

    theory = parse_tecplot(paths["u+y+theory.dat"])
    if theory["zones"][0]["data"].shape[1] < 5:
        raise RuntimeError("u+y+theory.dat has unexpected column count")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="CBS3D SA flat-plate paper plots with NASA/TMR references"
    )
    parser.add_argument(
        "--post-dir",
        required=True,
        help="post_flatplate_sa_global directory",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="output directory; default <post-dir>/paper_plots",
    )
    parser.add_argument(
        "--ref-dir",
        default=None,
        help="TMR reference cache; default <post-dir>/tmr_reference",
    )
    parser.add_argument(
        "--uinf",
        type=float,
        default=1.0,
        help="freestream velocity; default 1.0",
    )

    args = parser.parse_args()

    post_dir = os.path.abspath(args.post_dir)
    if not os.path.isdir(post_dir):
        raise RuntimeError("post directory not found: {}".format(post_dir))

    out_dir = (
        os.path.abspath(args.out_dir)
        if args.out_dir
        else os.path.join(post_dir, "paper_plots")
    )
    ref_dir = (
        os.path.abspath(args.ref_dir)
        if args.ref_dir
        else os.path.join(post_dir, "tmr_reference")
    )

    ensure_dir(out_dir)
    apply_paper_style()

    print("=" * 72)
    print("CBS3D SA FLAT-PLATE PUBLICATION POST-PROCESSING")
    print("=" * 72)
    print("post-dir : {}".format(post_dir))
    print("out-dir  : {}".format(out_dir))
    print("ref-dir  : {}".format(ref_dir))
    print("")

    paths = ensure_tmr_references(ref_dir)
    sanity_check_references(paths)
    print("[ref] official TMR reference parse: PASS")

    cbs_x, cbs_cf = load_cbs_cf(post_dir)
    profiles = load_cbs_profiles(post_dir)

    p1 = nearest_cbs_profile(profiles, TARGET_X1)
    p2 = nearest_cbs_profile(profiles, TARGET_X2)

    cbs_cf_x097 = interp_at(cbs_x, cbs_cf, TARGET_X1)

    refs_cf = ref_cf_distribution(paths)
    refs_conv = ref_cf_convergence(paths)

    plot_cf_distribution(cbs_x, cbs_cf, refs_cf, out_dir)
    plot_cf_convergence(cbs_cf_x097, refs_conv, out_dir)

    plot_uplus(p1, TARGET_X1, paths, out_dir)
    plot_uplus(p2, TARGET_X2, paths, out_dir)

    plot_velocity(p1, TARGET_X1, paths, out_dir, args.uinf)
    plot_velocity(p2, TARGET_X2, paths, out_dir, args.uinf)

    plot_mut(p1, paths, out_dir, args.uinf)
    plot_residual_history(post_dir, out_dir)

    write_validation_summary(
        out_dir,
        post_dir,
        cbs_cf_x097,
        profiles,
        paths,
    )

    print("")
    print("[CBS] x~0.97008 actual mesh station : {:.12g}".format(p1["x"]))
    print("[CBS] x~1.90334 actual mesh station : {:.12g}".format(p2["x"]))
    print("[CBS] first real y+ @ x~0.97008     : {:.12g}".format(first_real_yplus(p1)))
    print("[CBS] first real y+ @ x~1.90334     : {:.12g}".format(first_real_yplus(p2)))
    print("[CBS] Cf(x=0.97008)                 : {:.12g}".format(cbs_cf_x097))

    same = ref_incompressible_same_grid(paths)
    for title in sorted(same.keys()):
        ref = same[title]["cf"]
        err = percent_error(cbs_cf_x097, ref)
        print(
            "[TMR] {:35s} Cf={:.12g}  CBS error={:+.3f}%".format(
                title, ref, err
            )
        )

    print("")
    print("=" * 72)
    print("PAPER POST-PROCESSING: PASS")
    print("=" * 72)

    for path in sorted(glob.glob(os.path.join(out_dir, "*.png"))):
        print(os.path.basename(path))

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(2)
