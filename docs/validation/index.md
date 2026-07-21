# Validation programme

## Terminology

The project uses the following terms strictly:

- **Built** — the source compiled and linked.
- **Ran** — the executable advanced a case without a fatal runtime error.
- **Verified** — the implemented equations, discretisation or software behaviour were checked against a known numerical expectation.
- **Validated** — relevant output was compared with trusted analytical, experimental or published physical data.

A converged residual alone is not validation.

## Evidence classes

### 1. Legacy regression

The modern C++ solver is compared with the corresponding legacy Fortran implementation using matched:

```text
mesh
boundary conditions
material properties
initial conditions
time-step controls
solver tolerances
output sampling
```

This establishes behavioural continuity but does not independently validate the physics.

### 2. Independent benchmarks

At least one trusted reference is required for every released physical capability. Suitable references include analytical solutions, published benchmark data and experimental measurements.

### 3. Discretisation verification

The effect of numerical resolution must be quantified using:

- spatial mesh refinement;
- physical or pseudo-time-step refinement where relevant;
- pressure-solver tolerance sensitivity;
- iteration/convergence sensitivity;
- domain-size sensitivity where the boundary location affects the solution.

### 4. Parallel verification

The distributed solver must reproduce the accepted serial result before speed-up is reported. See [parallel status](../parallel/status.md).

## Validation matrix

| Capability | Primary case | Comparison quantity | Current status |
|---|---|---|---|
| Incompressible laminar flow | Lid-driven cavity, Re = 100 | Centreline velocity profiles | Existing development evidence; formal archive required |
| Internal viscous flow | Channel/Poiseuille flow | Analytical velocity profile and pressure gradient | Planned formal case |
| Transient vortex shedding | Flow around a cylinder | Strouhal number, drag/lift histories | Planned after transient path is accepted |
| Energy transport | Thermal diffusion/advection case | Analytical or manufactured temperature field | Planned formal case |
| Conjugate heat transfer | Fluid-solid benchmark | Interface temperature and energy balance | Development case available; formal reference required |
| Spalart–Allmaras | Turbulent flat plate | Skin friction, profiles and eddy viscosity | Active development |
| Distributed solver | Accepted serial cases | Field norms and conservation | Not yet accepted |

## Legacy benchmark lineage

Zeta Computational Resources provides public descriptions and downloadable benchmark data for the legacy CBSFlow family, including:

- 3D lid-driven cavity;
- 2D/3D flow around a cylinder;
- 2D lid-driven cavity;
- backward-facing step.

These data are useful for provenance and regression planning. Copyright and redistribution restrictions must be respected. Link to the original material rather than copying protected packages into this repository without written permission.

- [Zeta 3D benchmark data](https://www.zetacomp.com/benchmarks/benchmarks3d.asp)
- [Zeta 2D benchmark data](https://www.zetacomp.com/benchmarks/benchmarks2d.asp)

## Minimum case record

Every accepted validation result must record:

```text
case identifier
Git commit
input checksums
mesh statistics
material properties
boundary conditions
solver configuration
compiler and dependencies
hardware and parallel layout
convergence criteria
wall time
reference source
comparison method
error measures
plots and raw extracted values
```

## Required error measures

For a numerical field \(q\) compared at common points:

\[
E_2 = \left(
\frac{\sum_i(q_i-q_i^{\mathrm{ref}})^2}
{\sum_i(q_i^{\mathrm{ref}})^2}
\right)^{1/2},
\]

\[
E_\infty = \max_i |q_i-q_i^{\mathrm{ref}}|.
\]

Also report physically relevant integral quantities such as:

- net mass imbalance;
- drag and lift coefficients;
- pressure drop;
- wall shear or skin friction;
- integrated heat flux;
- fluid-solid energy imbalance.

## Acceptance practice

Acceptance tolerances must be specified before a new implementation is compared. A case should include:

1. machine-readable input;
2. documented reference data;
3. extraction script;
4. comparison script;
5. threshold definition;
6. archived summary table;
7. representative visualisation.

## Release gate

A capability can be marked **validated** in the root README only when its evidence is reproducible from repository or formally archived data, and the corresponding reference is cited.
