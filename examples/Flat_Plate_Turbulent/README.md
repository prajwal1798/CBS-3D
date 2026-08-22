# Turbulent zero-pressure-gradient flat plate

## Purpose

This directory is the canonical verification case for the Spalart-Allmaras
implementation and the production momentum wall treatment on
`integration/sa-mpi`.

The reference backbone is the NASA Langley Turbulence Modeling Resource (TMR)
2-D zero-pressure-gradient flat-plate SA case:

- plate starts at `x = 0` and extends to `x = 2`;
- upstream symmetry region begins at `x = -1/3`;
- `Re_x = 5.0e6` at `x = 1` and therefore `Re_x = 1.0e7` at `x = 2`;
- freestream SA working variable uses `nu_tilde_inf / nu_inf = 3`;
- the present CBS verification is incompressible and thin-span 3-D, with
  symmetry on both spanwise planes.

## Generate a solver-native case

No retained Sunbird input deck is required.  Generate `.plt`, `.bco`, `.par`
and a JSON mesh/condition manifest directly:

```bash
python3 examples/Flat_Plate_Turbulent/generate_case.py \
  --output-prefix flat_plate_yplus1 \
  --target-yplus 1
```

The generator uses a structured brick lattice and a globally conforming
six-TET4 Freudenthal split.  Wall-normal spacing is geometric.  Its first
height is estimated from

```text
Cf(x=1) = 0.0592 Re_x^(-1/5)
y1      = yplus_target * nu_inf / u_tau
u_inf   = 1 / Re_x
u_tilde_inf = 3 nu_inf
```

The reported target is a mesh-design value; accepted `y+` must be recomputed
from the converged wall traction.

## Production wall-treatment activation

The wall model is deliberately opt-in:

```bash
CBS3D_SA_WALL_TREATMENT=1 ./cbs3dpp_si flat_plate_yplus1
```

Optional smooth-wall constants are exposed only for verification studies:

```bash
CBS3D_WALL_KAPPA=0.41
CBS3D_WALL_B=5.2
```

With `CBS3D_SA_WALL_TREATMENT` unset, the solver retains the established
wall-resolved no-slip momentum boundary condition.

## What the production coupling changes

On fluid BC 530/532 faces only, the momentum wall treatment:

1. keeps impermeability as a strong constraint using the span of all incident
   wall normals;
2. leaves tangent-space velocity unconstrained;
3. removes the natural viscous Step-1 face contribution already assembled by
   the P1 operator;
4. replaces it by the constant-TRI3 Spalding traction
   `A_f/3 * (tau_w/rho)` at each face node;
5. rejects CHT activation until thermal/interface wall treatment is separately
   derived and verified.

The wall-normal projector is based on the local normal span, not an averaged
normal.  A planar wall removes one velocity DOF, an edge can remove two, and a
rank-three corner removes all three.

## Validation ladder

### A. Coupling regression: mandatory before any CFD comparison

CI must pass all of the following:

- Spalding algebra and wide-`Re_y` root solve;
- viscous-sublayer asymptote;
- TRI3 wall-load conservation;
- proof that the old natural viscous flux is *replaced*, not double counted;
- dissipative wall work;
- planar/edge/corner normal-span impermeability;
- full serial regression suite;
- full MPI/PETSc production build.

### B. Near-wall bridge case: `y+ ~ 1`

This is the first physical wall-model comparison because the current strong
SA boundary condition

```text
nu_tilde_wall = 0
```

is still consistent with a near-wall first sample.  Compare wall-model ON and
wall-resolved momentum on the same or matched meshes.  Required quantities:

- local `Cf(x)` over `0 < x <= 2`;
- integrated viscous drag;
- `y+(x)` and `u_tau(x)`;
- velocity profiles in inner variables;
- `nu_tilde`, `nu_t`, production, destruction and diffusion bounds;
- pressure uniformity along the plate;
- residual histories and mass conservation.

The wall-model bridge is accepted only if the skin-friction distribution and
integrated drag converge toward the wall-resolved SA reference without a
systematic traction sign or magnitude bias.

### C. Coarse-wall study: `y+ = 5, 10, 30, 50, 100`

**Do not treat these as model-validation results yet.**  The present
`nu_tilde=0` wall condition belongs to wall-resolved SA.  Moving the first
sample into the buffer/log layer requires a model-consistent SA scalar wall
boundary/treatment.  That is the next turbulence-model stage after the
`y+~1` momentum bridge.

Once that scalar treatment is implemented, regenerate the mesh family with:

```bash
for yp in 5 10 30 50 100; do
  python3 examples/Flat_Plate_Turbulent/generate_case.py \
    --output-prefix "flat_plate_yplus${yp}" \
    --target-yplus "${yp}"
done
```

For every level report `Cf`, drag, `y+`, velocity profiles, SA fields,
convergence and mesh sensitivity.  A successful coarse-wall model must show a
broad `y+` plateau rather than one accidentally favourable spacing.

## Reference comparison

Use NASA TMR SA flat-plate results as the primary numerical-verification
reference.  A simple turbulent flat-plate correlation may be plotted as a
secondary engineering check, but it is not a substitute for the TMR SA
solution or grid-convergence evidence.
