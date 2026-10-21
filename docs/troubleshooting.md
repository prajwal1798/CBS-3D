# Troubleshooting

Diagnose the earliest incorrect stage. Do not tune solver tolerances before confirming that the executable, case files and mesh interpretation are correct.

## Emergency checklist

1. Confirm the intended simulation type.
2. Confirm all required files exist and share one base name.
3. Confirm the executable supports the selected pressure solver.
4. Read the complete startup audit.
5. Check mesh, boundary and material totals.
6. Check the first non-finite or divergent quantity.
7. Preserve the log and input checksums before modifying the case.

## Build failures

### OpenMP unavailable

Symptoms:

```text
OpenMP not found
unknown option -fopenmp
```

Actions:

- verify that the compiler supports OpenMP;
- use a compiler-specific flag through CMake rather than adding flags manually;
- check that configuration did not silently disable OpenMP;
- clean the build directory after changing compilers.

### MPI wrapper missing

```text
mpicxx: command not found
```

Load or install the intended MPI implementation, then confirm:

```bash
which mpicxx
mpicxx --version
mpicxx --showme 2>/dev/null || mpicxx -show
```

### PETSc headers or libraries missing

Check:

```bash
printf '%s\n' "$PETSC_DIR"
find "$PETSC_DIR" -name petsc.h -o -name petscconf.h
find "$PETSC_DIR" -name 'libpetsc*'
```

`PETSC_DIR`, `PETSC_ARCH`, the compiler and MPI implementation must describe the same PETSc installation.

### Runtime library not found

```text
libpetsc.so: cannot open shared object file
```

Inspect:

```bash
ldd ./cbs3dpp_si_petsc
```

Correct the module environment, rpath or runtime library path. Do not copy arbitrary shared libraries beside the executable.

## Input failures

### Cannot open case file

The solver argument is the case base name.

For:

```text
input/cavity.plt
input/cavity.bco
input/cavity.par
```

run from the repository root with:

```bash
./cbs3dpp_si input/cavity
```

or from `input/` with:

```bash
../cbs3dpp_si cavity
```

### Invalid mesh header

Verify that the first `.plt` line contains the expected counts:

```text
NELEM NPOIN NBOUN
```

Then count connectivity, coordinate and boundary records.

### Boundary raw ID not found

Every raw boundary tag in `.plt` must be mapped in `.bco`:

```text
raw_boundary_id solver_boundary_id
```

Regenerate or correct the boundary mapping; do not change a physical tag blindly to bypass validation.

### Boundary face has no parent tetrahedron

The exporter failed to associate a boundary triangle with its adjacent tetrahedron. Rebuild the face-to-element map and verify orientation before running.

### Material header or connectivity mismatch

The `.material` element count and connectivity must match `.plt` exactly. Regenerate the material file from the same `.msh` input.

### Turbulence run requests material mode

Spalart–Allmaras requires molecular kinematic viscosity:

\[
\nu=\mu/\rho.
\]

Enable dimensional material properties and provide `.material` and `.matprop`.

## Pressure-solver failures

### Native CG reaches the iteration limit

Check in this order:

1. pressure reference or null-space treatment;
2. pressure boundary classification;
3. positive element volumes;
4. correct time-step/operator scaling;
5. finite pressure RHS;
6. diagonal/preconditioner values;
7. matrix-vector implementation;
8. tolerance and maximum iterations.

Do not increase the iteration limit before excluding a singular or incorrectly assembled system.

### PETSc option selected with the wrong binary

A PETSc pressure option requires a PETSc-enabled executable such as:

```text
cbs3dpp_si_petsc
```

The reference `cbs3dpp_si` binary cannot provide PETSc merely because runtime PETSc flags are appended.

### PETSc/HYPRE divergence

Run with diagnostic options:

```bash
-ksp_converged_reason
-ksp_monitor
-ksp_view
```

Record the operator size, KSP type, PC type, initial residual and divergence reason. Compare with the native operator on a small case.

## Non-physical solution

### Velocity is zero at the inlet

Check:

- boundary tag and `.bco` mapping;
- inlet area and mass-flow calculation;
- inlet value in `.par`;
- strong boundary application after Step 1 and Step 3;
- whether inlet nodes were also classified as wall/interface nodes.

### Pressure oscillates or diverges

Check time-step size, pressure reference, mesh quality, boundary consistency, reversed tetrahedra and velocity boundary leakage.

### Temperature rises on the wrong surface

Check:

- raw and mapped heat-flux IDs;
- parent boundary element and outward normal;
- flux units and sign convention;
- fluid-solid material labels;
- face area calculation;
- whether the intended face is visually selected in the boundary VTU.

### SA variable remains zero

Check:

- turbulence activation and model selection;
- positive molecular viscosity;
- inlet/farfield SA value;
- wall SA condition;
- wall-distance field;
- production, diffusion and destruction diagnostics;
- clipping or non-finite-value handling.

### Boundary layer begins before the plate leading edge

The upstream lower boundary was likely classified as no-slip. Use symmetry/slip upstream and no-slip only from the physical leading edge onward.

## Parallel failures

### Global totals change with rank count

Possible causes:

- duplicate ownership;
- missing ghost entities;
- incorrect reduction;
- label loss during distribution;
- counting local ghosts as global owned entities.

### Rank-local values disagree after assembly

Check the exact communication sequence:

```text
forward input halo
owned-element assembly
reverse additive halo
owner update
boundary enforcement
forward updated halo
```

### MPI build runs but only rank zero advances

This is the current bootstrap behaviour, not evidence of distributed acceleration. Use the phase acceptance tests in [parallel status](parallel/status.md).

## Required bug report data

Include:

```text
Git commit
platform and compiler
build command
executable name
case-file checksums
complete startup output
first failing iteration
pressure solver and options
mesh/material/boundary totals
stack trace or sanitizer output
smallest reproducible case
```
