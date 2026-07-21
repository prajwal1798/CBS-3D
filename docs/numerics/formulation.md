# Numerical formulation

## Scope

The active CBS3D++ path is a semi-implicit incompressible pressure-projection solver on linear tetrahedral finite elements. It is a practical CBS specialisation rather than a literal implementation of every term in the most general Split-A formulation.

The accepted current scope is:

- incompressible flow;
- constant-property or element-wise constant material properties;
- steady/pseudo-time iteration with fixed or locally evaluated time-step controls;
- optional energy transport over fluid and solid regions;
- linear P1 tetrahedral elements;
- lumped nodal mass;
- native matrix-free CG or PETSc pressure solution.

True physical-time transient integration is not yet a released capability.

## Governing equations

For an incompressible Newtonian fluid:

\[
\nabla\cdot\mathbf{u}=0,
\]

\[
\rho\left(
\frac{\partial\mathbf{u}}{\partial t}
+\mathbf{u}\cdot\nabla\mathbf{u}
\right)
=-\nabla p+\nabla\cdot\boldsymbol{\tau}+\rho\mathbf{b}.
\]

The active laminar momentum diffusion is implemented as a component-wise Laplacian form. It should not be described as the complete deviatoric-stress operator unless that formulation is explicitly assembled.

For the thermal field:

\[
\rho c_p\left(
\frac{\partial T}{\partial t}
+\mathbf{u}\cdot\nabla T
\right)
=\nabla\cdot(k\nabla T)+Q.
\]

In solid elements, velocity is zero and the active equation reduces to heat diffusion with the relevant material properties.

## P1 tetrahedral discretisation

Each tetrahedron has four linear shape functions. Their spatial gradients are constant within the element:

\[
\nabla N_a = \text{constant in } \Omega_e.
\]

The element volume follows from the coordinate Jacobian:

\[
V_e=\frac{|\det J_e|}{6}.
\]

The implementation stores one-based scientific arrays to preserve the legacy indexing model. Index zero is intentionally unused in the principal solver containers.

## CBS update sequence

### Step 1: momentum predictor

A pressure-free intermediate velocity is formed from convection, diffusion and active source terms:

\[
\mathbf{u}^{*}=\mathbf{u}^{n}+\Delta\mathbf{u}^{*}.
\]

Element contributions are accumulated to nodes and divided by the active lumped momentum mass/time diagonal. Strong velocity boundary conditions are reapplied after the nodal update.

### Step 2: pressure equation

The pressure-Poisson system enforces the incompressibility constraint. In abstract form:

\[
\mathbf{K}_p\,p^{n+1}=\mathbf{r}_p(\mathbf{u}^{*}).
\]

The native path applies the pressure operator matrix-free using element stiffness data. The PETSc path assembles the corresponding algebraic system for KSP solution.

A pressure reference or equivalent null-space treatment is required for purely Neumann pressure configurations.

### Step 3: velocity correction

The pressure gradient corrects the predicted velocity:

\[
\mathbf{u}^{n+1}
=\mathbf{u}^{*}
-\Delta t\,\rho^{-1}\nabla p^{n+1}.
\]

The active implementation omits the higher-order characteristic pressure contribution present in some general CBS derivations. This is a documented formulation choice and must be preserved or deliberately revised with regression evidence.

### Step 4: energy/scalar transport

When enabled, the temperature field is advanced after the corrected fluid velocity:

\[
T^{n+1}=T^n+\Delta T.
\]

Fluid elements include advection and diffusion. Solid elements include diffusion and active sources/flux terms. Conformal fluid-solid interface nodes provide temperature continuity through the shared nodal field.

## Mass matrices

The solver uses lumped nodal mass terms for explicit or semi-explicit field updates. The pressure operator is assembled from the tetrahedral gradient products:

\[
K_{ab}^{(e)}
=V_e\,\nabla N_a\cdot\nabla N_b,
\]

with the active time and density scaling applied by the pressure-system preparation path.

## Boundary conditions

Boundary preprocessing maps raw physical tags to solver boundary identifiers, determines parent-element local faces, computes outward normals and classifies boundary nodes.

Strong boundary values are reapplied at the mathematically appropriate points in the four-stage update. Changes to boundary sequencing can alter the projection and must be treated as numerical changes, not cosmetic refactoring.

## Pressure solvers

### Native

- matrix-free pressure multiplication;
- preconditioned conjugate gradient;
- optional Jacobi-style diagonal preconditioning;
- residual and iteration reporting.

### PETSc

- KSP pressure solution;
- PETSc runtime options;
- HYPRE/BoomerAMG when available and requested.

The present PETSc production path is not yet a distributed pressure solve over rank-owned degrees of freedom.

## Turbulence

Spalart–Allmaras development introduces a transported modified turbulent viscosity and eddy-viscosity contribution. It requires dimensional molecular viscosity and density. Its implementation status and validation evidence must be reported separately from the laminar solver.

## Numerical-change policy

Any change to the following requires a documented regression:

- sign or scaling of pressure RHS/operator terms;
- time-step placement;
- mass lumping;
- boundary application order;
- shape gradients or Jacobian orientation;
- material-domain masks;
- pressure null-space treatment;
- fluid-solid interface assembly;
- turbulence source, diffusion or destruction terms.

The detailed file-by-file audit remains under `docs/CBS3Dpp_SI_Full_Solver_Audit.md` and should be updated when the active formulation changes.
