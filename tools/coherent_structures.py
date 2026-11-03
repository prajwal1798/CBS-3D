#!/usr/bin/env pvpython
# =============================================================================
# CBS3D++_SI
#
# FILE      : tools/coherent_structures.py
# PURPOSE   : PUBLICATION FIGURES FROM ONE CBS3D SNAPSHOT (DISTRIBUTED .pvd /
#             .pvtu OUTPUT) WITH PARAVIEW 6.1 (pvpython, NO GUI NEEDED).
#
#   --what structures   VORTEX / COHERENT-STRUCTURE IDENTIFICATION (DEFAULT)
#   --what fields       TEMPERATURE AND VELOCITY FIELD FIGURES
#   --what all          BOTH
#
#   STRUCTURES, FOR EACH VELOCITY FIELD ("velocity" AND, IF PRESENT, THE
#   TIME-AVERAGED "mean_velocity"), ON THE FLUID ELEMENTS ONLY:
#
#     J        = grad u,   J_ij = du_i/dx_j                (VTK Gradient filter)
#     S, W     = SYMMETRIC / ANTISYMMETRIC PARTS OF J
#     omega    = curl u,   |omega| VORTICITY MAGNITUDE
#     Q        = 0.5 ( |W|^2 - |S|^2 ),  Q > 0 IN VORTEX CORES
#                (HUNT, WRAY & MOIN 1988, CTR SUMMER PROGRAM CTR-S88)
#     lambda2  = MIDDLE EIGENVALUE OF S^2 + W^2,  < 0 IN VORTEX CORES
#                (JEONG & HUSSAIN 1995, J FLUID MECH 285:69-94)
#     lambdaci = |IMAG PART| OF THE COMPLEX EIGENVALUE PAIR OF J ("SWIRLING
#                STRENGTH"),  > 0 WHERE THE FLUID SWIRLS
#                (ZHOU, ADRIAN, BALACHANDAR & KENDALL 1999, J FLUID MECH
#                 387:353-396)
#     Hn       = u.omega / (|u| |omega|) IN [-1, 1], NORMALISED HELICITY. ITS
#                SIGN GIVES THE SENSE OF ROTATION ABOUT THE LOCAL FLOW
#                DIRECTION, SO THE TWO CELLS OF A DEAN-VORTEX PAIR BEHIND A
#                BEND APPEAR IN OPPOSITE COLOURS
#                (LEVY, DEGANI & SEGINER 1990, AIAA J 28:1347-1352)
#
#   ISO-LEVELS ARE SET RELATIVE TO A ROBUST REFERENCE, THE 99.5th PERCENTILE
#   OF THE POSITIVE VALUES IN THE FLUID, SO ONE FRACTION f WORKS FOR ALL
#   MESHES AND TIMES:
#     Q_iso = f Q_ref,   lambda2_iso = -f l2_ref,
#     lambdaci_iso = sqrt(f) lci_ref     (lambdaci SCALES LIKE 1/s, Q LIKE 1/s^2)
#
# MERGING   : THE PIECES OF A DISTRIBUTED SNAPSHOT SHARE THEIR INTERFACE NODES
#             (vtkGhostType = 1 ON NON-OWNED COPIES). THEY ARE MERGED FIRST
#             (CLEAN TO GRID). WITHOUT THIS THE GRADIENT IS ONE-SIDED AT EVERY
#             PARTITION BOUNDARY (SPURIOUS SPECKLES IN Q) AND TRANSLUCENT
#             SURFACES SHOW THE PARTITION FACES AS A NETWORK OF LINES.
#
# OUTPUT    : <out>/*.png                   FIGURES (WHITE BACKGROUND)
#             <out>/coherent_structures.pvsm   PARAVIEW STATE, FILE > LOAD STATE
#             <out>/summary.json            REFERENCE VALUES AND ISO-LEVELS
#             <out>/fluid_criteria.vtu      (--save-vtu) FLUID FIELD WITH ALL
#                                           CRITERIA, FOR LATER ANALYSIS
#
# USAGE     : pvpython --force-offscreen-rendering coherent_structures.py \
#                 ~/PhD/CFD/runs/run30082/blanket_30082.pvd
#             OPTIONS: --what structures|fields|all   --out DIR
#                      --fraction 0.05   --sweep 0.01 0.02 0.05 0.1 0.2
#                      --fields velocity mean_velocity   --fluid-material 0
#                      --solid-material 1   --velocity-label "time-averaged velocity"
#                      --resolution 2400 1350   --no-clean   --save-vtu
#                      --plain-text (NO LaTeX-STYLE LABELS)   --raytrace
# =============================================================================

import argparse
import builtins  # pvpython PRELOADS paraview.simple *, WHICH SHADOWS sum/min/max
import json
import math
import os
import sys

from paraview.simple import (
    CleantoGrid, ColorBy, Contour, CreateView, ExtractSurface, GetColorTransferFunction,
    GetOpacityTransferFunction, GetScalarBar, Gradient, Hide, OpenDataFile, PassArrays,
    ProgrammableFilter, Render, SaveData, SaveScreenshot, SaveState, Show, Slice, Text, Threshold)


# =============================================================================
# PROGRAMMABLE-FILTER BODY: lambda2, lambdaci, HELICITY, CONSISTENCY CHECKS
# =============================================================================
CRITERIA_SCRIPT = r'''
import numpy as np
FIELDS = %(fields)r
NAMES = %(names)r
inp = inputs[0]
for F in FIELDS:
    J = np.asarray(inp.PointData[F + "_grad"], dtype=float).reshape(-1, 3, 3)
    om_vtk = np.asarray(inp.PointData[F + "_vorticity"], dtype=float)
    # ================================
    # COMPONENT ORDER CHECK: WITH J_ij = du_i/dx_j THE CURL IS
    # (J21 - J12, J02 - J20, J10 - J01). IF THE TRANSPOSE WAS STORED THE SIGN
    # FLIPS; THE RESULT IS THEN INDEPENDENT OF THE STORAGE CONVENTION.
    # ================================
    om_J = np.stack([J[:, 2, 1] - J[:, 1, 2], J[:, 0, 2] - J[:, 2, 0], J[:, 1, 0] - J[:, 0, 1]], axis=1)
    scale = np.abs(om_vtk).max() + 1e-300
    if np.abs(om_J - om_vtk).max() > 1e-6 * scale:
        if np.abs(om_J + om_vtk).max() <= 1e-6 * scale:
            J = np.transpose(J, (0, 2, 1))
        else:
            raise RuntimeError("velocity gradient and VTK vorticity disagree for " + F)
    Jt = np.transpose(J, (0, 2, 1))
    S = 0.5 * (J + Jt)
    W = 0.5 * (J - Jt)
    Q = 0.5 * ((W * W).sum(axis=(1, 2)) - (S * S).sum(axis=(1, 2)))
    l2 = np.linalg.eigvalsh(S @ S + W @ W)[:, 1]
    lci = np.abs(np.linalg.eigvals(J).imag).max(axis=1)
    U = np.asarray(inp.PointData[F], dtype=float)
    om_mag = np.sqrt((om_vtk * om_vtk).sum(axis=1))
    u_mag = np.sqrt((U * U).sum(axis=1))
    hel = (U * om_vtk).sum(axis=1)
    small = (u_mag < 1e-12) | (om_mag < 1e-12)
    Hn = np.where(small, 0.0, hel / (u_mag * om_mag + 1e-300))
    q_vtk = np.asarray(inp.PointData[F + "_Q"], dtype=float)
    q_err = np.abs(Q - q_vtk).max() / (np.abs(q_vtk).max() + 1e-300)
    output.PointData.append(om_mag, F + "_vorticity_magnitude")
    output.PointData.append(l2, F + "_lambda2")
    output.PointData.append(lci, F + "_lambdaci")
    output.PointData.append(hel, F + "_helicity")
    output.PointData.append(Hn, F + "_normalised_helicity")
    output.PointData.append(u_mag, F + "_magnitude")
    def p995(a):
        return float(np.percentile(a, 99.5)) if a.size else 0.0
    ref = [p995(Q[Q > 0]), p995(-l2[l2 < 0]), p995(lci[lci > 0]),
           float(np.percentile(om_mag, 1.0)), float(np.percentile(om_mag, 99.9)),
           float(Q.max()), float(-l2.min()), float(lci.max()), float(q_err), float(u_mag.max())]
    for name, value in zip(NAMES, ref):
        output.FieldData.append(np.array([value]), F + "_ref_" + name)
'''

REFERENCE_NAMES = ["Q_ref", "lambda2_ref", "lambdaci_ref", "vorticity_p01", "vorticity_p999",
                   "Q_max", "minus_lambda2_max", "lambdaci_max", "Q_relative_check", "speed_max"]


# =============================================================================
# LABELS (MATHTEXT NEEDS MATPLOTLIB, WHICH THE KITWARE BINARIES INCLUDE)
# =============================================================================
class Labels(object):
    def __init__(self, plain):
        self.plain = plain

    def __call__(self, math_text, plain_text):
        return plain_text if self.plain else math_text


# =============================================================================
# SMALL HELPERS
# =============================================================================
def field_data_value(proxy, name):
    """READ A ONE-VALUE FIELD-DATA ARRAY FROM THE DATA INFORMATION (NO DATA COPY)."""
    info = proxy.GetDataInformation().DataInformation.GetFieldDataInformation()
    for i in range(info.GetNumberOfArrays()):
        ai = info.GetArrayInformation(i)
        if ai.GetName() == name:
            return ai.GetComponentRange(0)[0]
    raise RuntimeError("field data array %s not found" % name)


def point_array_names(proxy):
    pinfo = proxy.GetDataInformation().DataInformation.GetPointDataInformation()
    return [pinfo.GetArrayInformation(i).GetName() for i in range(pinfo.GetNumberOfArrays())]


def point_array_range(proxy, name):
    pinfo = proxy.GetDataInformation().DataInformation.GetPointDataInformation()
    for i in range(pinfo.GetNumberOfArrays()):
        ai = pinfo.GetArrayInformation(i)
        if ai.GetName() == name:
            return ai.GetComponentRange(-1 if ai.GetNumberOfComponents() > 1 else 0)
    return None


def bounds(proxy):
    return proxy.GetDataInformation().GetBounds()


def set_camera(view, b, mode):
    """'iso': OBLIQUE PERSPECTIVE VIEW, THINNEST AXIS TOWARDS THE VIEWER.
       'top': PARALLEL PROJECTION ALONG THE THINNEST AXIS (NO TRAPEZOID DISTORTION)."""
    ext = [b[1] - b[0], b[3] - b[2], b[5] - b[4]]
    c = [0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3]), 0.5 * (b[4] + b[5])]
    n = ext.index(builtins.min(ext))
    a, bb = [i for i in range(3) if i != n]
    up = [0.0, 0.0, 0.0]
    d = [0.0, 0.0, 0.0]
    up[bb] = 1.0
    if mode == "top":
        d[n] = 1.0
    else:
        d[a], d[bb], d[n] = 0.8, 0.35, 1.0
    r = math.sqrt(builtins.sum(e * e for e in ext))
    view.CameraFocalPoint = c
    view.CameraPosition = [c[i] + 3.0 * r * d[i] for i in range(3)]
    view.CameraViewUp = up
    view.CameraParallelProjection = 1 if mode == "top" else 0
    view.ResetCamera(False)
    if mode == "top":
        # OBJECT FILLS ~78 % OF THE WIDTH, LEAVING ROOM FOR THE COLOUR BAR
        aspect = float(view.ViewSize[0]) / float(view.ViewSize[1])
        view.CameraParallelScale = 0.5 * builtins.max(ext[a] / (0.78 * aspect), ext[bb] / 0.85)
        # MOVE THE CAMERA ALONG SCREEN-RIGHT SO THE OBJECT SITS LEFT OF THE BAR
        look = [-d[i] for i in range(3)]
        right = [look[1] * up[2] - look[2] * up[1], look[2] * up[0] - look[0] * up[2],
                 look[0] * up[1] - look[1] * up[0]]
        shift = [0.08 * ext[a] * right[i] for i in range(3)]
        view.CameraFocalPoint = [c[i] + shift[i] for i in range(3)]
        view.CameraPosition = [c[i] + shift[i] + 3.0 * r * d[i] for i in range(3)]
    return n, c


def style_view(view, raytrace):
    view.OrientationAxesVisibility = 1
    view.OrientationAxesLabelColor = [0.1, 0.1, 0.1]
    view.UseColorPaletteForBackground = 0
    view.Background = [1.0, 1.0, 1.0]
    try:
        view.BackgroundColorMode = "Single Color"
    except AttributeError:
        pass
    if raytrace:
        try:
            view.EnableRayTracing = 1
            view.BackEnd = "OSPRay pathtracer"
            view.SamplesPerPixel = 16
            view.Denoise = 1
        except AttributeError:
            print("  NOTE: ray tracing not available in this ParaView build; using OpenGL")


def apply_preset(lut, name):
    """PRESET NAMES CHANGED IN PARAVIEW 6 ('Viridis (matplotlib)' -> 'Viridis')."""
    for n in (name, name + " (matplotlib)"):
        try:
            lut.ApplyPreset(n, True)
            return
        except RuntimeError:
            pass
    print("  NOTE: colour map %s not found; using the default" % name)


def colour(display, array, cmap, crange, log=False):
    ColorBy(display, ("POINTS", array))
    lut = GetColorTransferFunction(array)
    apply_preset(lut, cmap)
    lut.AutomaticRescaleRangeMode = "Never"
    lut.UseLogScale = 0
    lut.RescaleTransferFunction(crange[0], crange[1])
    if log:
        lut.MapControlPointsToLogSpace()
        lut.UseLogScale = 1
    GetOpacityTransferFunction(array).RescaleTransferFunction(crange[0], crange[1])
    return lut


def scalar_bar(view, lut, title, fmt="{:.4g}"):
    bar = GetScalarBar(lut, view)
    bar.Title = title
    bar.ComponentTitle = ""
    bar.TitleColor = [0.1, 0.1, 0.1]
    bar.LabelColor = [0.1, 0.1, 0.1]
    bar.TitleFontSize = 20
    bar.LabelFontSize = 16
    bar.ScalarBarLength = 0.4
    bar.WindowLocation = "Any Location"
    bar.Position = [0.88, 0.08]
    for prop, val in (("RangeLabelFormat", fmt), ("LabelFormat", fmt), ("AutomaticLabelFormat", 0)):
        try:
            setattr(bar, prop, val)
        except AttributeError:
            pass
    bar.Visibility = 1
    return bar


def title_text(view, text):
    t = Text(Text=text)
    d = Show(t, view)
    d.FontSize = 20
    d.Color = [0.1, 0.1, 0.1]
    d.WindowLocation = "Upper Left Corner"
    return t


class Figures(object):
    """ONE RENDER VIEW; EACH FIGURE SHOWS ITS OWN PROXIES AND HIDES THEM AGAIN."""

    def __init__(self, view, out, res):
        self.view, self.out, self.res = view, out, res
        self.shown, self.bars = [], []

    def show(self, proxy):
        d = Show(proxy, self.view)
        self.shown.append(proxy)
        return d

    def save(self, name):
        Render(self.view)
        SaveScreenshot(os.path.join(self.out, name), self.view, ImageResolution=self.res)
        for p in self.shown:
            Hide(p, self.view)
        for b in self.bars:
            b.Visibility = 0
        self.shown, self.bars = [], []
        print("   wrote", name)


# =============================================================================
# MAIN
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description="CBS3D publication figures (ParaView 6.1)")
    ap.add_argument("input", help=".pvd or .pvtu of ONE snapshot")
    ap.add_argument("--what", choices=["structures", "fields", "all"], default="structures")
    ap.add_argument("--out", default=None, help="output folder (default: <input folder>/figures)")
    ap.add_argument("--fields", nargs="+", default=["velocity", "mean_velocity"])
    ap.add_argument("--fluid-material", type=int, default=0)
    ap.add_argument("--solid-material", type=int, default=1)
    ap.add_argument("--velocity-label", default=None,
                    help='description of the "velocity" array, e.g. "time-averaged velocity" for an energy run')
    ap.add_argument("--fraction", type=float, default=0.05, help="iso-level as fraction of the reference")
    ap.add_argument("--sweep", type=float, nargs="*", default=[0.01, 0.02, 0.05, 0.1, 0.2],
                    help="extra Q fractions rendered to choose a level (give no values to skip)")
    ap.add_argument("--resolution", type=int, nargs=2, default=[2400, 1350])
    ap.add_argument("--no-clean", action="store_true", help="do not merge the partition pieces")
    ap.add_argument("--save-vtu", action="store_true")
    ap.add_argument("--plain-text", action="store_true", help="plain ASCII labels instead of MathText")
    ap.add_argument("--raytrace", action="store_true")
    args = ap.parse_args()

    L = Labels(args.plain_text)
    src_path = os.path.abspath(os.path.expanduser(args.input))
    out = os.path.abspath(os.path.expanduser(args.out)) if args.out else \
        os.path.join(os.path.dirname(src_path), "figures")
    os.makedirs(out, exist_ok=True)
    res = list(args.resolution)
    do_struct = args.what in ("structures", "all")
    do_fields = args.what in ("fields", "all")

    # ================================
    # READ, KEEP THE NEEDED ARRAYS, MERGE PIECES
    # ================================
    print("== reading", src_path)
    reader = OpenDataFile(src_path)
    if reader is None:
        sys.exit("cannot open " + src_path)
    reader.UpdatePipeline()
    names = point_array_names(reader)
    fields = [f for f in args.fields if f in names]
    if not fields:
        sys.exit("none of %s found; point arrays are %s" % (args.fields, names))
    has_T = "temperature" in names
    info = reader.GetDataInformation()
    print("   points %d  cells %d  velocity fields %s  temperature %s"
          % (info.GetNumberOfPoints(), info.GetNumberOfCells(), fields, "yes" if has_T else "no"))

    keep = fields + [n for n in ("temperature",) if n in names]
    slim = PassArrays(Input=reader, PointDataArrays=keep, CellDataArrays=["material_id"])
    merged = slim if args.no_clean else CleantoGrid(Input=slim)

    def material(lo, hi):
        return Threshold(Input=merged, Scalars=["CELLS", "material_id"],
                         LowerThreshold=lo, UpperThreshold=hi, ThresholdMethod="Between")

    fluid = material(args.fluid_material, args.fluid_material)
    fluid.UpdatePipeline()
    finfo = fluid.GetDataInformation()
    print("   fluid points %d  cells %d" % (finfo.GetNumberOfPoints(), finfo.GetNumberOfCells()))

    lab = {"velocity": args.velocity_label or "instantaneous velocity",
           "mean_velocity": "time-averaged velocity"}
    summary = {"input": src_path, "fraction": args.fraction, "fields": {}}

    view = CreateView("RenderView")
    view.ViewSize = res
    style_view(view, args.raytrace)
    fig = Figures(view, out, res)

    # ================================
    # ONE CAMERA FRAME FOR EVERY FIGURE: THE WHOLE SLICE (STEEL + COOLANT).
    # MID-PLANES ARE NORMAL TO ITS THINNEST AXIS, THROUGH THE COOLANT CENTRE.
    # ================================
    domain = material(builtins.min(args.fluid_material, args.solid_material),
                      builtins.max(args.fluid_material, args.solid_material))
    domain.UpdatePipeline()
    b_dom = bounds(domain)
    b_flu = bounds(fluid)
    ext = [b_dom[1] - b_dom[0], b_dom[3] - b_dom[2], b_dom[5] - b_dom[4]]
    n_axis = ext.index(builtins.min(ext))
    normal = [1.0 if i == n_axis else 0.0 for i in range(3)]
    centre = [0.5 * (b_flu[0] + b_flu[1]), 0.5 * (b_flu[2] + b_flu[3]), 0.5 * (b_flu[4] + b_flu[5])]

    def plane(inp):
        s_ = Slice(Input=inp)
        s_.SliceType = "Plane"
        s_.SliceType.Origin = centre
        s_.SliceType.Normal = normal
        return s_

    walls = ExtractSurface(Input=fluid)

    # ================================
    # COHERENT STRUCTURES
    # ================================
    crit = None
    if do_struct:
        last = fluid
        for f in fields:
            last = Gradient(Input=last, ScalarArray=["POINTS", f], ComputeGradient=1,
                            ComputeVorticity=1, ComputeQCriterion=1, ResultArrayName=f + "_grad",
                            VorticityArrayName=f + "_vorticity", QCriterionArrayName=f + "_Q")
        crit = ProgrammableFilter(Input=last, CopyArrays=1,
                                  Script=CRITERIA_SCRIPT % {"fields": fields, "names": REFERENCE_NAMES})
        crit.UpdatePipeline()

        for f in fields:
            ref = dict((n, field_data_value(crit, f + "_ref_" + n)) for n in REFERENCE_NAMES)
            if ref["Q_relative_check"] > 1e-6:
                print("   WARNING: Q from VTK and from the gradient differ by %.2e (relative)"
                      % ref["Q_relative_check"])
            ref["Q_iso"] = args.fraction * ref["Q_ref"]
            ref["lambda2_iso"] = -args.fraction * ref["lambda2_ref"]
            ref["lambdaci_iso"] = math.sqrt(args.fraction) * ref["lambdaci_ref"]
            summary["fields"][f] = ref
            print("== %s: Q_ref %.4g 1/s^2 (max %.4g) | -lambda2_ref %.4g | lambdaci_ref %.4g 1/s | check %.1e"
                  % (f, ref["Q_ref"], ref["Q_max"], ref["lambda2_ref"], ref["lambdaci_ref"],
                     ref["Q_relative_check"]))
            print("   iso-levels at f = %g: Q = %.4g, lambda2 = %.4g, lambdaci = %.4g"
                  % (args.fraction, ref["Q_iso"], ref["lambda2_iso"], ref["lambdaci_iso"]))

        if args.save_vtu:
            path = os.path.join(out, "fluid_criteria.vtu")
            SaveData(path, proxy=crit, CompressorType="ZLib")
            print("   wrote", path)

        mid = plane(crit)

        speed = L(r"$\vert\mathbf{u}\vert$  [m s$^{-1}$]", "|u| [m/s]")

        def iso(f, array, level, col, cmap, crange, bar_title, name, title):
            cont = Contour(Input=crit, ContourBy=["POINTS", array], Isosurfaces=[level])
            d = fig.show(cont)
            lut = colour(d, col, cmap, crange)
            d.Specular = 0.3
            d.SpecularPower = 40.0
            wd = fig.show(walls)
            ColorBy(wd, None)
            wd.AmbientColor = wd.DiffuseColor = [0.6, 0.62, 0.66]
            wd.Opacity = 0.12
            fig.bars.append(scalar_bar(view, lut, bar_title))
            fig.shown.append(title_text(view, title))
            set_camera(view, b_dom, "iso")
            fig.save(name)

        def midplane(array, cmap, crange, log, bar_title, name, title):
            d = fig.show(mid)
            lut = colour(d, array, cmap, crange, log)
            fig.bars.append(scalar_bar(view, lut, bar_title))
            fig.shown.append(title_text(view, title))
            set_camera(view, b_dom, "top")
            fig.save(name)

        for f in fields:
            ref = summary["fields"][f]
            vmax = ref["speed_max"] or 1.0
            lb = lab.get(f, f)
            if ref["Q_ref"] > 0:
                for s in [args.fraction] + [x for x in args.sweep if x != args.fraction]:
                    name = f + "_Q.png" if s == args.fraction else "%s_Q_sweep_%g.png" % (f, s)
                    iso(f, f + "_Q", s * ref["Q_ref"], f + "_magnitude", "Viridis", [0.0, vmax], speed, name,
                        L(r"$Q$-criterion, %s:  $Q = %.3g\ \mathrm{s^{-2}}$  ($%.3g\,Q_\mathrm{ref}$)"
                          % (lb, s * ref["Q_ref"], s),
                          "Q-criterion, %s: Q = %.3g 1/s^2 (%.3g Q_ref)" % (lb, s * ref["Q_ref"], s)))
            if ref["lambda2_ref"] > 0:
                iso(f, f + "_lambda2", ref["lambda2_iso"], f + "_magnitude", "Viridis", [0.0, vmax], speed,
                    f + "_lambda2.png",
                    L(r"$\lambda_2$-criterion, %s:  $\lambda_2 = %.3g\ \mathrm{s^{-2}}$" % (lb, ref["lambda2_iso"]),
                      "lambda2-criterion, %s: lambda2 = %.3g 1/s^2" % (lb, ref["lambda2_iso"])))
            if ref["lambdaci_ref"] > 0:
                iso(f, f + "_lambdaci", ref["lambdaci_iso"], f + "_normalised_helicity", "Cool to Warm",
                    [-1.0, 1.0], L(r"$H_n$", "normalised helicity"), f + "_lambdaci_helicity.png",
                    L(r"Swirling strength, %s:  $\lambda_{ci} = %.3g\ \mathrm{s^{-1}}$, coloured by "
                      r"normalised helicity $H_n$ (sense of rotation)" % (lb, ref["lambdaci_iso"]),
                      "Swirling strength, %s: lambda_ci = %.3g 1/s, coloured by normalised helicity"
                      % (lb, ref["lambdaci_iso"])))
            hi = builtins.max(ref["vorticity_p999"], 1e-12)
            lo = builtins.max(ref["vorticity_p01"], 1e-3 * hi)
            midplane(f + "_vorticity_magnitude", "Inferno", [lo, hi], True,
                     L(r"$\vert\boldsymbol{\omega}\vert$  [s$^{-1}$]", "|vorticity| [1/s]"),
                     f + "_vorticity_midplane.png",
                     L(r"Vorticity magnitude on the mid-plane, %s" % lb,
                       "Vorticity magnitude on the mid-plane, %s" % lb))
            midplane(f + "_normalised_helicity", "Cool to Warm", [-1.0, 1.0], False,
                     L(r"$H_n$", "normalised helicity"), f + "_helicity_midplane.png",
                     L(r"Normalised helicity $H_n$ on the mid-plane, %s" % lb,
                       "Normalised helicity on the mid-plane, %s" % lb))

    # ================================
    # TEMPERATURE AND VELOCITY FIELDS
    # ================================
    if do_fields:
        solid = material(args.solid_material, args.solid_material)

        def field_fig(proxy, array, cmap, crange, bar_title, name, title, mode, bb):
            d = fig.show(proxy)
            d.Opacity = 1.0
            lut = colour(d, array, cmap, crange)
            d.Specular = 0.2
            fig.bars.append(scalar_bar(view, lut, bar_title))
            fig.shown.append(title_text(view, title))
            set_camera(view, bb, mode)
            fig.save(name)

        if has_T:
            tr = point_array_range(domain, "temperature")
            tbar = L(r"$T$  [K]", "T [K]")
            s_dom, s_sol, s_flu = ExtractSurface(Input=domain), ExtractSurface(Input=solid), walls
            field_fig(s_dom, "temperature", "Inferno", tr, tbar, "temperature_iso.png",
                      "Temperature, whole slice (steel and coolant surfaces)", "iso", b_dom)
            field_fig(s_sol, "temperature", "Inferno", tr, tbar, "temperature_steel_iso.png",
                      "Temperature, steel only", "iso", b_dom)
            field_fig(s_flu, "temperature", "Inferno", tr, tbar, "temperature_coolant_iso.png",
                      "Temperature, coolant only (channel walls)", "iso", b_dom)
            field_fig(plane(domain), "temperature", "Inferno", tr, tbar, "temperature_midplane_top.png",
                      "Temperature on the mid-plane", "top", b_dom)
            field_fig(s_dom, "temperature", "Inferno", tr, tbar, "temperature_top.png",
                      "Temperature, outer surface, viewed normal to the slice", "top", b_dom)
        for f in fields:
            vr = point_array_range(fluid, f)
            vmax = vr[1] if vr else 1.0
            d_mid = plane(fluid)
            field_fig(d_mid, f, "Viridis", [0.0, vmax], L(r"$\vert\mathbf{u}\vert$  [m s$^{-1}$]", "|u| [m/s]"),
                      f + "_magnitude_midplane_top.png",
                      "Speed on the coolant mid-plane, %s" % lab.get(f, f), "top", b_dom)

    # ================================
    # LEAVE A USEFUL STATE FOR THE GUI
    # ================================
    if crit is not None:
        f0 = fields[0]
        if summary["fields"][f0]["Q_ref"] > 0:
            cont = Contour(Input=crit, ContourBy=["POINTS", f0 + "_Q"],
                           Isosurfaces=[summary["fields"][f0]["Q_iso"]])
            d = Show(cont, view)
            colour(d, f0 + "_magnitude", "Viridis", [0.0, summary["fields"][f0]["speed_max"] or 1.0])
            wd = Show(walls, view)
            ColorBy(wd, None)
            wd.Opacity = 0.12
        set_camera(view, b_dom, "iso")
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    state = os.path.join(out, "coherent_structures.pvsm")
    SaveState(state)
    print("== done:", out)
    print("   GUI: File > Load State >", state)


if __name__ == "__main__":
    main()
