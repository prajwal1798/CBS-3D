# CBS3D++_SI Spalart-Allmaras Implementation Plan

Status: design-freeze draft 0.1  
Scope: serial/OpenMP first, MPI-compatible data model second, PETSc pressure solve retained separately  
Primary reference: NASA Turbulence Modeling Resource, Spalart-Allmaras model page: https://tmbwg.github.io/turbmodels/spalart.html

---

## 1. Purpose

This document defines the disciplined implementation plan for adding the Spalart-Allmaras one-equation turbulence model to CBS3D++_SI.

The objective is not merely to add an eddy-viscosity variable. The objective is to implement a verifiable RANS extension of the existing semi-implicit incompressible CBS finite-element solver while preserving the validated laminar solver path.

The turbulence implementation must satisfy these requirements:

1. Laminar mode must remain numerically unchanged when turbulence is disabled.
2. The turbulence model must be explicitly selectable and documented.
3. Wall distance must be computed as true minimum distance to physical no-slip wall faces, not by grid-line search or nearest wall node.
4. Turbulent viscosity must be applied only in fluid momentum equations.
5. Turbulent thermal diffusivity must be applied only in fluid energy equations.
6. Solid heat conduction must not receive turbulent conductivity.
7. Fluid-solid CHT interface velocity/turbulence behaviour must be treated as a wall for turbulence, while temperature continuity remains through the thermal system.
8. PETSc must remain responsible only for the pressure linear solve unless a later scalar-transport backend is deliberately introduced.
9. The first benchmark must be a canonical verification case before any blanket or fusion geometry is trusted.

---

## 2. Model version selected

### 2.1 Baseline model

The baseline model is the NASA TMR standard Spalart-Allmaras model, denoted here as `SA`.

The model solves the transported working variable:

```text
nu_tilde
```

The turbulent dynamic viscosity is then computed from:

```text
mu_t = rho * nu_tilde * f_v1
```

For incompressible constant-density regions this is equivalently:

```text
nu_t = nu_tilde * f_v1
mu_t = rho * nu_t
```

where:

```text
nu = mu / rho
chi = nu_tilde / nu
```

### 2.2 Robustness option

The implementation must be structured so that `SA_NEG` can be enabled later without rewriting the solver.

The default first validation should be the non-negative standard `SA` path. However, the code should not be designed in a way that prevents the later SA-neg branch.

Configuration names:

```text
turbulence_model = SA
turbulence_model = SA_NEG
```

Initial implementation target:

```text
SA only, with robust limiting and positivity protection
```

Later extension target:

```text
SA_NEG negative branch
```

---

## 3. Governing equations

### 3.1 Existing incompressible momentum equation

CBS3D++_SI currently advances the incompressible velocity field using the CBS fractional-step method. In laminar mode, the viscous stress uses the molecular dynamic viscosity:

```text
mu = mu_e
```

With SA enabled, the effective fluid viscosity becomes:

```text
mu_eff = mu + mu_t
```

where `mu_t` is computed from the SA model.

The turbulent viscosity must not be applied in solid elements.

### 3.2 SA transport equation

The standard SA model is written in non-conservative form as:

```text
D(nu_tilde)/Dt = Production - Destruction + Diffusion + NonlinearDiffusion
```

In index notation:

```text
partial(nu_tilde)/partial(t)
+ u_j partial(nu_tilde)/partial(x_j)
=
cb1 * (1 - ft2) * S_tilde * nu_tilde
-
[cw1 * fw - cb1 / kappa^2 * ft2] * (nu_tilde / d)^2
+
1 / sigma * {
    partial/partial(x_j)[(nu + nu_tilde) partial(nu_tilde)/partial(x_j)]
    + cb2 * partial(nu_tilde)/partial(x_i) * partial(nu_tilde)/partial(x_i)
}
```

where:

```text
d = minimum distance to nearest physical no-slip wall
```

### 3.3 Eddy viscosity relation

```text
nu_t = nu_tilde * fv1
mu_t = rho * nu_t
```

with:

```text
fv1 = chi^3 / (chi^3 + cv1^3)
chi = nu_tilde / nu
```

### 3.4 Auxiliary functions

The standard SA auxiliary definitions are:

```text
S_tilde = Omega + nu_tilde / (kappa^2 d^2) * fv2
```

```text
fv2 = 1 - chi / (1 + chi * fv1)
```

```text
fw = g * [(1 + cw3^6) / (g^6 + cw3^6)]^(1/6)
```

```text
g = r + cw2 * (r^6 - r)
```

```text
r = min[nu_tilde / (S_tilde * kappa^2 * d^2), 10]
```

```text
ft2 = ct3 * exp(-ct4 * chi^2)
```

The vorticity magnitude is:

```text
Omega = sqrt(2 * W_ij * W_ij)
```

where:

```text
W_ij = 0.5 * (du_i/dx_j - du_j/dx_i)
```

### 3.5 S_tilde limiting

The calculation of `S_tilde` must be protected against zero or negative values before computing `r`.

The implementation shall use the Allmaras-Johnson-Spalart style limiter documented by NASA TMR:

Define:

```text
S_bar = nu_tilde / (kappa^2 d^2) * fv2
```

Then:

```text
if S_bar >= -c2 * Omega:
    S_tilde = Omega + S_bar
else:
    S_tilde = Omega + Omega * (c2^2 * Omega + c3 * S_bar) / ((c3 - 2*c2) * Omega - S_bar)
```

with:

```text
c2 = 0.7
c3 = 0.9
```

If `S_tilde == 0` during the computation of `r`, set:

```text
r = 10
```

This limiting choice must be printed in the run summary.

---

## 4. Model constants

Use the constants from the NASA TMR standard SA formulation:

| Constant | Value |
|---|---:|
| `cb1` | 0.1355 |
| `sigma` | 2.0 / 3.0 |
| `cb2` | 0.622 |
| `kappa` | 0.41 |
| `cw2` | 0.3 |
| `cw3` | 2.0 |
| `cv1` | 7.1 |
| `ct3` | 1.2 |
| `ct4` | 0.5 |

The remaining constant is derived:

```text
cw1 = cb1 / kappa^2 + (1 + cb2) / sigma
```

Thermal turbulent Prandtl number for fluid heat transfer:

```text
Pr_t = 0.90
```

This value is used only when turbulent thermal diffusivity is enabled.

---

## 5. Boundary conditions

### 5.1 No-slip physical wall

For no-slip physical walls:

```text
nu_tilde = 0
nu_t = 0
mu_t = 0
```

This applies to:

```text
solid-fluid wall/interface for momentum and turbulence
external no-slip walls
```

### 5.2 Inlet

For fully turbulent inlet/farfield behaviour:

```text
nu_tilde_inlet = sa_inlet_ratio * nu_molecular
```

Initial accepted values:

```text
sa_inlet_ratio = 3.0 to 5.0
```

Default:

```text
sa_inlet_ratio = 3.0
```

The selected value must be printed in the run summary.

### 5.3 Outlet

Use zero normal gradient:

```text
d(nu_tilde)/dn = 0
```

For the first implementation, this means no direct overwrite at outlet nodes after the update unless a later weak boundary term is explicitly implemented.

### 5.4 Solid region

No SA equation is solved in solid elements.

At nodes that are solid-only:

```text
nu_tilde = 0
nu_t = 0
mu_t = 0
```

At fluid-solid interface nodes:

```text
SA boundary condition: wall value nu_tilde = 0
thermal behaviour: shared temperature continuity remains unchanged
```

### 5.5 Partition boundaries

Partition boundaries are artificial MPI communication interfaces.

They are not physical walls and must not receive SA wall boundary conditions.

---

## 6. Wall distance

Wall distance is a required field:

```text
wall_distance(node) = d(node)
```

The value must be the true minimum Euclidean distance from the node to the nearest physical no-slip wall triangle.

Forbidden approximations:

```text
nearest wall node distance
searching along grid lines
assuming wall-normal grid lines
using partition boundary faces
```

Serial first implementation:

```text
for every active fluid node:
    d = min distance to all physical wall triangles
```

Production implementation:

```text
AABB tree or bounding-volume hierarchy over wall triangles
```

Validation test:

```text
flat wall at y = 0:
    expected wall_distance = y
```

---

## 7. Discretisation strategy

### 7.1 First implementation

The first implementation shall use the existing CBS finite-element infrastructure:

```text
P1 tetrahedral basis
lumped nodal update
local element assembly
existing element geometry dNkdx and detJ
existing OpenMP colouring where appropriate
```

### 7.2 SA update form

The first stable implementation is a segregated scalar transport update:

```text
nu_tilde_new = nu_tilde_old + dt * M_lumped^{-1} * RHS_SA
```

where `RHS_SA` includes:

```text
advection
diffusion
nonlinear gradient term
production
destruction
```

For stability, the first version may treat stiff destruction semi-implicitly by placing a positive destruction coefficient in the denominator:

```text
nu_tilde_new = (nu_tilde_old + dt * explicit_positive_terms) / (1 + dt * destruction_coefficient)
```

The exact split must be documented in code comments and in this document after implementation.

### 7.3 Positivity protection

For the initial `SA` implementation:

```text
nu_tilde = max(nu_tilde, 0)
```

This is a numerical protection for the non-negative model.

When `SA_NEG` is implemented later, negative values are handled by the SA-neg branch and turbulent eddy viscosity is set to zero for negative `nu_tilde`.

---

## 8. Coupling order inside one solver iteration

Recommended segregated order:

```text
1. Copy old velocity, pressure, temperature and nu_tilde fields.
2. Compute nu_t, mu_t, mu_eff and k_eff from current nu_tilde.
3. CBS Step 1: momentum predictor using mu_eff in fluid elements.
4. CBS Step 2: pressure solve using PETSc KSPCG path when enabled.
5. CBS Step 3: velocity correction.
6. SA scalar transport solve using corrected velocity.
7. Recompute nu_t, mu_t, mu_eff and k_eff.
8. CBS Step 4: energy solve using k_eff in fluid elements only.
9. Compute laminar and turbulence residuals.
10. Write output fields.
```

The SA field is therefore lagged with respect to the momentum solve by one nonlinear iteration, which is acceptable for the first segregated RANS implementation.

---

## 9. Effective fluid properties

### 9.1 Momentum viscosity

For fluid elements:

```text
mu_eff_e = mu_e + mu_t_e
```

For solid elements:

```text
mu_eff_e = mu_e
```

Momentum assembly must use:

```text
mu_used = turbulence_on && fluid_element ? mu_eff_e : mu_e
```

### 9.2 Thermal conductivity

For fluid elements with turbulent heat transfer enabled:

```text
k_eff_e = k_e + rho_e * cp_e * nu_t_e / Pr_t
```

For solid elements:

```text
k_eff_e = k_e
```

Energy assembly must use:

```text
k_used = turbulence_on && fluid_element ? k_eff_e : k_e
```

This prevents turbulent conductivity from contaminating solid heat conduction.

---

## 10. Required new state variables

Add to `CBSStateSI`:

```cpp
Array1D<Real> nu_tilde;
Array1D<Real> nu_tilde1;
Array1D<Real> nu_t;
Array1D<Real> mu_t;
Array1D<Real> wall_distance;
Array1D<Real> sa_rhs;
Array1D<Real> sa_residual;
Array1D<Real> sa_source;
Array1D<Real> sa_mass;

Array1D<Real> nu_tilde_e;
Array1D<Real> nu_t_e;
Array1D<Real> mu_t_e;
Array1D<Real> mu_eff_e;
Array1D<Real> k_eff_e;
```

Do not replace:

```cpp
Array1D<Real> rho_e;
Array1D<Real> cp_e;
Array1D<Real> k_e;
Array1D<Real> mu_e;
```

Those remain molecular/material properties.

---

## 11. Required configuration additions

Add to `RunConfig`:

```cpp
Int turbulence_on = 0;
Int turbulence_model = 0;          // 0 none, 1 SA, 2 SA_NEG
Real sa_inlet_ratio = 3.0;
Real sa_prandtl_t = 0.90;
Real sa_min_wall_distance = 1.0e-14;
Int sa_use_stilde_limiter = 1;
Int sa_clip_negative = 1;
Int turbulent_heat_transfer = 1;
Int sa_write_diagnostics = 1;
```

The parser must accept text aliases:

```text
NONE
SA
SA_NEG
```

---

## 12. Required new source modules

Add:

```text
include/cbs/turbulence/SpalartAllmaras.hpp
src/turbulence/SpalartAllmaras.cpp

include/cbs/turbulence/WallDistance.hpp
src/turbulence/WallDistance.cpp

include/cbs/assembly/SpalartAllmarasAssembly.hpp
src/assembly/SpalartAllmarasAssembly.cpp

include/cbs/boundary/TurbulenceBoundary.hpp
src/boundary/TurbulenceBoundary.cpp
```

Later optional module:

```text
include/cbs/io/TurbulenceInput.hpp
src/io/TurbulenceInput.cpp
```

---

## 13. Existing files to modify

```text
CMakeLists.txt
include/cbs/core/CBSStateSI.hpp
include/cbs/core/RunConfig.hpp
src/io/MeshIO.cpp
src/preprocess/Preprocess.cpp
src/boundary/Boundary.cpp
src/timestep/TimeStep.cpp
src/assembly/MomentumAssembly.cpp
src/assembly/EnergyAssembly.cpp
src/solver/Steps.cpp
src/solver/Solver.cpp
src/io/Post.cpp
```

---

## 14. PETSc scope

PETSc remains responsible for:

```text
pressure KSPCG solve
pressure preconditioner
pressure matrix/vector storage when PETSc is enabled
```

PETSc is not used initially for:

```text
SA scalar transport
wall-distance computation
eddy-viscosity update
momentum assembly
energy assembly
boundary conditions
```

This keeps turbulence implementation transparent and debuggable.

---

## 15. MPI compatibility plan

SA is first implemented in serial/OpenMP mode.

The MPI data model must later treat SA fields like scalar nodal fields.

Owner-to-ghost update:

```text
nu_tilde
nu_t
mu_t
wall_distance if needed
```

Ghost-to-owner accumulation:

```text
sa_rhs
sa_mass
sa_residual sum
```

Element-level turbulent properties remain local:

```text
nu_tilde_e
nu_t_e
mu_t_e
mu_eff_e
k_eff_e
```

Partition boundaries must not become turbulence walls.

---

## 16. Output additions

VTU output should include:

```text
nu_tilde
nu_t
mu_t
mu_eff
k_eff
wall_distance
chi
SA residual if enabled
```

Console summary should include:

```text
turbulence model
SA inlet ratio
SA Pr_t
minimum wall distance
maximum nu_t / nu
maximum mu_t / mu
SA residual
```

---

## 17. Verification and validation ladder

### 17.1 Unit tests

Test pure SA functions:

```text
chi
fv1
fv2
S_tilde limiter
r limiter
fw
ft2
production term
destruction term
eddy viscosity conversion
```

### 17.2 Wall-distance verification

Flat-wall geometry:

```text
wall at y = 0
expected d = y
```

### 17.3 Scalar transport manufactured test

Before full RANS coupling, verify:

```text
advection
diffusion
source assembly
boundary condition application
```

using a known artificial field.

### 17.4 NASA TMR flat plate verification

Primary benchmark:

```text
2D zero-pressure-gradient flat plate with SA
```

Quantities to compare:

```text
skin-friction coefficient Cf(x)
eddy viscosity profiles
maximum eddy viscosity versus x
velocity profile in wall units
integrated drag coefficient if applicable
```

NASA TMR provides grids and expected SA reference results for this case.

Important limitation:

```text
NASA reference case is compressible M = 0.2.
CBS3D++_SI is incompressible.
Agreement should be close in low-Mach attached boundary-layer quantities, but not bitwise or exactly identical.
```

### 17.5 Turbulent pipe/channel validation

After flat-plate verification, validate an internal-flow case:

```text
fully developed turbulent pipe
or turbulent channel
```

Compare:

```text
friction factor
wall shear stress
pressure gradient
mean velocity profile
u+ vs y+
```

---

## 18. No-break rule for laminar CBS3D++_SI

When:

```text
turbulence_on = 0
```

then:

```text
mu_eff_e = mu_e
k_eff_e = k_e
nu_tilde arrays are inactive
SA boundary conditions are skipped
SA residual is not part of convergence
existing laminar validation must remain unchanged
```

Any commit that changes laminar lid-driven cavity, pipe or blanket results without an explicit reason is rejected.

---

## 19. Implementation milestones

### Milestone 1: documentation freeze

This file is created and reviewed before code begins.

### Milestone 2: state and configuration

Add fields and parser controls. Compile only. No numerical effect.

### Milestone 3: wall distance

Implement true point-to-triangle wall distance. Verify on flat wall.

### Milestone 4: SA algebra

Implement pure functions and unit tests.

### Milestone 5: SA boundary conditions

Apply wall, inlet, outlet and solid-region SA states.

### Milestone 6: SA transport assembly

Implement scalar transport update for `nu_tilde`.

### Milestone 7: eddy-viscosity coupling

Use `mu_eff_e` in momentum and `k_eff_e` in fluid energy.

### Milestone 8: OpenMP safety

Make SA assembly race-free with element colouring or a verified alternative.

### Milestone 9: flat-plate verification

Compare against NASA TMR SA reference data.

### Milestone 10: pipe/channel validation

Validate turbulent internal flow.

### Milestone 11: MPI compatibility

Add SA halo update and reverse accumulation after serial/OpenMP correctness.

---

## 20. Immediate next action

No code implementation should begin until this document is accepted as the active design basis.

After acceptance, the next commit should perform only:

```text
RunConfig additions
CBSStateSI field additions
CMakeLists additions for empty turbulence modules
no change to laminar numerical output
```
