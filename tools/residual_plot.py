#!/usr/bin/env python3
"""
Publication figures: CBS3D++_SI pressure-Poisson solver comparison.

Compares the two pressure-Poisson equation (PPE) solvers on the same case:
    PCG(Jacobi) : matrix-free conjugate gradient, Jacobi (diagonal) preconditioner
    PCG(AMG)    : PETSc KSPCG, algebraic-multigrid preconditioner
                  (BoomerAMG via hypre, or PETSc GAMG - check with -ksp_view)

Terminology used on all axes (standard CFD / numerical linear algebra):
    time step n        : one CBS pseudo-time step of the steady-state march
    normalised residual: relative change-based convergence residual R
    PPE iterations     : Krylov iterations of the pressure-Poisson solve

Figures produced (PDF + PNG, captionless - captions belong in the paper):
    fig_ppe_convergence_history   overlaid normalised-residual histories
                                  (demonstrates identical outer convergence:
                                   the preconditioner changes cost, not the answer)
    fig_ppe_iterations_per_step   PPE Krylov iterations per time step (semilogy)
                                  (the robustness/conditioning result)
    fig_ppe_cumulative_work       cumulative PPE iterations vs time step
                                  (total algorithmic work to convergence)
    table_ppe_solver_summary.csv  numerical summary + work-reduction factor

Expected input files inside --residual-dir:
    <case>_residuals_jacobian.csv   (PCG-Jacobi run)
    <case>_residuals_petsc.csv      (PCG-AMG run)

Usage:
    python residual_plot.py --residual-dir <dir> [--case LidDrivenCavity3D]
                            [--tolerance 1e-5] [--out-dir figs]
"""

from __future__ import annotations

from pathlib import Path
import argparse
import math
import shutil

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, LogLocator, NullFormatter


RESIDUAL_COLUMNS = ["u_rel", "v_rel", "w_rel", "p_rel"]

LABEL_JACOBI = r"PCG(Jacobi)"
LABEL_AMG = r"PCG(AMG)"


def configure_matplotlib(use_tex: bool = False) -> None:
    """Compact serif figures sized for two-column journal layout."""
    latex_available = shutil.which("latex") is not None
    plt.rcParams.update({
        "text.usetex": bool(use_tex and latex_available),
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "CMU Serif",
                       "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "figure.dpi": 180,
        "savefig.dpi": 500,
        "savefig.bbox": "tight",
        "axes.grid": True,
        "grid.alpha": 0.20,
        "grid.linewidth": 0.35,
        "axes.linewidth": 0.75,
        "lines.linewidth": 0.95,
        "legend.frameon": False,
    })


def read_residual_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required CSV: {path}")

    df = pd.read_csv(path)
    required = {
        "iteration", "u_rel", "v_rel", "w_rel", "p_rel",
        "cg_iterations", "cg_initial_l2", "cg_final_l2", "cg_relative_l2",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")

    df = df.copy()
    df["iteration"] = df["iteration"].astype(int)
    df = df.sort_values("iteration").drop_duplicates("iteration", keep="last")
    df = df.reset_index(drop=True)
    df["Rmax"] = df[RESIDUAL_COLUMNS].max(axis=1)
    df["cum_ppe_iterations"] = df["cg_iterations"].cumsum()
    return df


def trim_to_first_convergence(
    df: pd.DataFrame, tolerance: float
) -> tuple[pd.DataFrame, int | None]:
    mask = (df[RESIDUAL_COLUMNS] <= tolerance).all(axis=1)
    if not mask.any():
        return df.copy(), None

    stop_iteration = int(df.loc[mask, "iteration"].iloc[0])
    trimmed = df[df["iteration"] <= stop_iteration].copy().reset_index(drop=True)
    trimmed["cum_ppe_iterations"] = trimmed["cg_iterations"].cumsum()
    return trimmed, stop_iteration


def sample_for_clean_line(df: pd.DataFrame, every: int) -> pd.DataFrame:
    """First row, every Nth row, final row - line cleanliness only."""
    if df.empty or every <= 1:
        return df.copy()
    keep = (df["iteration"] == df["iteration"].iloc[0]) | \
           (df["iteration"] % every == 0) | \
           (df["iteration"] == df["iteration"].iloc[-1])
    return df.loc[keep].copy().reset_index(drop=True)


def clean_axes(ax) -> None:
    ax.tick_params(direction="out", length=3.0, width=0.7)
    for spine in ax.spines.values():
        spine.set_linewidth(0.75)


def savefig(fig: plt.Figure, outbase: Path) -> None:
    fig.savefig(outbase.with_suffix(".pdf"))
    fig.savefig(outbase.with_suffix(".png"))
    plt.close(fig)


# ============================================================================
# Figure 1: outer convergence history (both solvers overlap -> same physics)
# ============================================================================
def plot_component_convergence(
    df_full: pd.DataFrame,
    solver_label: str,
    outbase: Path,
    tolerance: float,
    sample_every: int,
) -> None:
    """Per-variable normalised residual history: R_u, R_v, R_w, R_p."""
    df = sample_for_clean_line(df_full, sample_every)

    series = [
        ("u_rel", r"$R_u$", "-"),
        ("v_rel", r"$R_v$", "--"),
        ("w_rel", r"$R_w$", "-."),
        ("p_rel", r"$R_p$", ":"),
    ]

    fig, ax = plt.subplots(figsize=(3.45, 2.30))
    for column, label, style in series:
        ax.semilogy(df["iteration"], df[column], linestyle=style, label=label)
    ax.axhline(tolerance, linestyle=":", linewidth=0.85, color="0.35")

    ax.set_xlabel(r"Time step $n$")
    ax.set_ylabel(r"Normalised residual")

    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_minor_formatter(NullFormatter())
    clean_axes(ax)
    ax.legend(loc="upper right", ncol=2, title=solver_label, title_fontsize=7)
    savefig(fig, outbase)


def plot_convergence_history(
    jacobi_full: pd.DataFrame,
    amg_full: pd.DataFrame,
    outbase: Path,
    tolerance: float,
    sample_every: int,
) -> None:
    jacobi = sample_for_clean_line(jacobi_full, sample_every)
    amg = sample_for_clean_line(amg_full, sample_every)

    fig, ax = plt.subplots(figsize=(3.45, 2.30))
    ax.semilogy(jacobi["iteration"], jacobi["Rmax"], label=LABEL_JACOBI)
    ax.semilogy(amg["iteration"], amg["Rmax"], linestyle="--", label=LABEL_AMG)
    ax.axhline(tolerance, linestyle=":", linewidth=0.85, color="0.35",
               label=rf"$R = 10^{{{int(math.log10(tolerance))}}}$")

    ax.set_xlabel(r"Time step $n$")
    ax.set_ylabel(r"Normalised residual $R$")

    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_minor_formatter(NullFormatter())
    clean_axes(ax)
    ax.legend(loc="upper right")
    savefig(fig, outbase)


# ============================================================================
# Figure 2: PPE Krylov iterations per time step (the conditioning result)
# ============================================================================
def plot_ppe_iterations_per_step(
    jacobi_full: pd.DataFrame,
    amg_full: pd.DataFrame,
    outbase: Path,
    sample_every: int,
) -> None:
    jacobi = sample_for_clean_line(jacobi_full, sample_every)
    amg = sample_for_clean_line(amg_full, sample_every)

    fig, ax = plt.subplots(figsize=(3.45, 2.30))
    ax.semilogy(jacobi["iteration"], jacobi["cg_iterations"], label=LABEL_JACOBI)
    ax.semilogy(amg["iteration"], amg["cg_iterations"], linestyle="--",
                label=LABEL_AMG)

    ax.set_xlabel(r"Time step $n$")
    ax.set_ylabel(r"PPE iterations per time step")

    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_minor_formatter(NullFormatter())
    clean_axes(ax)
    ax.legend(loc="center right")
    savefig(fig, outbase)


# ============================================================================
# Figure 3: cumulative PPE work to convergence
# ============================================================================
def plot_cumulative_work(
    jacobi_full: pd.DataFrame,
    amg_full: pd.DataFrame,
    outbase: Path,
    sample_every: int,
) -> None:
    jacobi = sample_for_clean_line(jacobi_full, sample_every)
    amg = sample_for_clean_line(amg_full, sample_every)

    fig, ax = plt.subplots(figsize=(3.45, 2.30))
    ax.semilogy(jacobi["iteration"], jacobi["cum_ppe_iterations"],
                label=LABEL_JACOBI)
    ax.semilogy(amg["iteration"], amg["cum_ppe_iterations"], linestyle="--",
                label=LABEL_AMG)

    ax.set_xlabel(r"Time step $n$")
    ax.set_ylabel(r"Cumulative PPE iterations")

    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_minor_formatter(NullFormatter())
    clean_axes(ax)
    ax.legend(loc="lower right")
    savefig(fig, outbase)


# ============================================================================
# Numerical summary table
# ============================================================================
def write_summary(
    jacobi: pd.DataFrame,
    amg: pd.DataFrame,
    outpath: Path,
    tolerance: float,
    jacobi_stop: int | None,
    amg_stop: int | None,
) -> None:
    rows = []
    for solver, df, stop in [
        (LABEL_JACOBI, jacobi, jacobi_stop),
        (LABEL_AMG, amg, amg_stop),
    ]:
        rows.append({
            "ppe_solver": solver,
            "convergence_tolerance": tolerance,
            "time_steps_to_convergence": stop if stop is not None else "not_reached",
            "total_ppe_iterations": int(df["cg_iterations"].sum()),
            "mean_ppe_iterations_per_step": float(df["cg_iterations"].mean()),
            "median_ppe_iterations_per_step": float(df["cg_iterations"].median()),
            "max_ppe_iterations_per_step": int(df["cg_iterations"].max()),
            "final_Ru": float(df["u_rel"].iloc[-1]),
            "final_Rv": float(df["v_rel"].iloc[-1]),
            "final_Rw": float(df["w_rel"].iloc[-1]),
            "final_Rp": float(df["p_rel"].iloc[-1]),
        })

    summary = pd.DataFrame(rows)
    j_total = summary.loc[0, "total_ppe_iterations"]
    a_total = summary.loc[1, "total_ppe_iterations"]
    j_mean = summary.loc[0, "mean_ppe_iterations_per_step"]
    a_mean = summary.loc[1, "mean_ppe_iterations_per_step"]
    summary["work_reduction_factor_total"] = [
        j_total / a_total if a_total else math.nan, 1.0]
    summary["work_reduction_factor_mean"] = [
        j_mean / a_mean if a_mean else math.nan, 1.0]
    summary.to_csv(outpath, index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--residual-dir", type=Path, required=True,
                        help="Directory containing the residual CSV files.")
    parser.add_argument("--case", type=str, default="LidDrivenCavity3D",
                        help="Case base name used in the CSV filenames.")
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--tolerance", type=float, default=1.0e-5)
    parser.add_argument("--sample-every", type=int, default=100,
                        help="Plot every Nth time step for clean lines.")
    parser.add_argument("--usetex", action="store_true")
    args = parser.parse_args()

    out_dir = args.out_dir if args.out_dir is not None \
        else args.residual_dir / "figures_paper"
    out_dir.mkdir(parents=True, exist_ok=True)

    configure_matplotlib(use_tex=args.usetex)

    jacobi_raw = read_residual_csv(
        args.residual_dir / f"{args.case}_residuals_jacobian.csv")
    amg_raw = read_residual_csv(
        args.residual_dir / f"{args.case}_residuals_petsc.csv")

    jacobi, jacobi_stop = trim_to_first_convergence(jacobi_raw, args.tolerance)
    amg, amg_stop = trim_to_first_convergence(amg_raw, args.tolerance)

    plot_component_convergence(
        jacobi, LABEL_JACOBI, out_dir / "fig_convergence_uvwp_jacobi",
        args.tolerance, args.sample_every)
    plot_component_convergence(
        amg, LABEL_AMG, out_dir / "fig_convergence_uvwp_amg",
        args.tolerance, args.sample_every)
    plot_convergence_history(
        jacobi, amg, out_dir / "fig_ppe_convergence_history",
        args.tolerance, args.sample_every)
    plot_ppe_iterations_per_step(
        jacobi, amg, out_dir / "fig_ppe_iterations_per_step",
        args.sample_every)
    plot_cumulative_work(
        jacobi, amg, out_dir / "fig_ppe_cumulative_work",
        args.sample_every)
    write_summary(
        jacobi, amg, out_dir / "table_ppe_solver_summary.csv",
        args.tolerance, jacobi_stop, amg_stop)

    j_tot = int(jacobi["cg_iterations"].sum())
    a_tot = int(amg["cg_iterations"].sum())
    print(f"{LABEL_JACOBI}: converged at time step {jacobi_stop}, "
          f"total PPE iterations {j_tot:,}")
    print(f"{LABEL_AMG}:    converged at time step {amg_stop}, "
          f"total PPE iterations {a_tot:,}")
    if a_tot:
        print(f"PPE work reduction: {j_tot / a_tot:.1f}x")
    print(f"Figures -> {out_dir}")


if __name__ == "__main__":
    main()
