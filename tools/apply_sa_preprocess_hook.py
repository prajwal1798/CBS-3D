from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT / "src" / "solver" / "Solver.cpp"

text = SOLVER.read_text(encoding="utf-8")

include_line = '#include "cbs/turbulence/TurbulencePreprocess.hpp"\n'
include_anchor = '#include "cbs/timestep/TimeStep.hpp"\n'

if include_line not in text:
    if include_anchor not in text:
        raise SystemExit("Solver.cpp include anchor not found")
    text = text.replace(include_anchor, include_anchor + include_line, 1)
    print("[patch] added TurbulencePreprocess include")
else:
    print("[skip] TurbulencePreprocess include already present")

old_line = "//     8. Apply the initial velocity, pressure and temperature conditions."
new_lines = """//     8. Apply the initial velocity, pressure and temperature conditions.
//     9. If requested, preprocess the Spalart-Allmaras turbulence model."""

if new_lines not in text:
    if old_line not in text:
        raise SystemExit("Solver.cpp initialise comment anchor not found")
    text = text.replace(old_line, new_lines, 1)
    print("[patch] updated initialise comment")
else:
    print("[skip] initialise comment already updated")

hook = """

        // Precompute all geometry-dependent Spalart-Allmaras quantities before
        // the CBS time loop.  In the current milestone this includes the
        // OpenMP wall-distance search and the initial eddy-viscosity field.
        if (s_.cfg.turbulence_on > 0)
        {
            Post::printStage(
                "Turbulence preprocessing",
                "Spalart-Allmaras wall distance");

            TurbulencePreprocess::prepareSpalartAllmaras(s_);

            Post::printStageDone(
                "Turbulence preprocessing",
                "SA wall distance ready");
        }
"""

anchor = '        Post::printStageDone("Initial boundary values", "initial field constrained");\n'

if "TurbulencePreprocess::prepareSpalartAllmaras" not in text:
    if anchor not in text:
        raise SystemExit("Solver.cpp initial-value anchor not found")
    text = text.replace(anchor, anchor + hook, 1)
    print("[patch] inserted SA preprocessing call")
else:
    print("[skip] SA preprocessing call already present")

SOLVER.write_text(text, encoding="utf-8")
print("SA preprocessing hook applied.")
