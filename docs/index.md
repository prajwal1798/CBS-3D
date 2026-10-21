<div align="center">

<img src="assets/cbs3d-logo.svg" alt="CBS3D++ — A Finite Element based Fluid Dynamics Solver" width="100%">

</div>

CBS3D++ is a C++20 finite-element solver for three-dimensional incompressible flow, multiphysics transport and high-performance CFD development. The active numerical method is a semi-implicit Characteristic-Based Split pressure-projection scheme on linear tetrahedral elements.

## Documentation map

| Area | Purpose |
|---|---|
| [Installation](getting-started/installation.md) | Build variants, dependencies and toolchain compatibility |
| [Quick start](getting-started/quickstart.md) | Minimal route from source to a first run |
| [Case files](user-guide/case-files.md) | Mesh, boundary, parameter and material-file contract |
| [Numerical formulation](numerics/formulation.md) | Governing equations, CBS sequence and discrete assumptions |
| [Parallel status](parallel/status.md) | Implemented infrastructure, missing distributed operations and acceptance gates |
| [Validation programme](validation/index.md) | Regression, benchmark and verification evidence |
| [Legacy regression](validation/legacy-regression.md) | Controlled comparison against the legacy Fortran solver |
| [Legacy provenance](legacy/provenance.md) | Numerical lineage, copyright boundary and redevelopment scope |
| [Sunbird workflow](hpc/sunbird.md) | Reproducible HPC build and execution policy |
| [Developer architecture](developer-guide/architecture.md) | Source layout, solver call order and change discipline |
| [Windows/MSYS2](developer-guide/windows-msys2.md) | Portable local development workflow |
| [Troubleshooting](troubleshooting.md) | Failure diagnosis by build stage and solver symptom |

## Solver status

The current accepted execution path is serial or OpenMP element assembly with either the native pressure solver or PETSc KSP/HYPRE. MPI mesh-distribution and halo-exchange components are present, but the complete CBS update is not yet distributed over ranks.

The documentation uses three status terms:

- **Available** — implemented in the active solver path and exercised by current cases.
- **Development** — implemented partially or undergoing numerical validation.
- **Planned** — design target without a completed accepted implementation.

Capabilities must not be advertised as released until the corresponding validation gate is recorded.

## Scientific scope

The current development programme covers:

- incompressible laminar flow;
- semi-implicit pressure projection;
- conjugate heat transfer and scalar energy transport;
- Spalart–Allmaras turbulence-model development;
- PETSc pressure-solver integration;
- OpenMP shared-memory acceleration;
- distributed ownership, halo exchange and domain decomposition;
- ParaView output and residual monitoring.

## Documentation policy

The root `README.md` is intentionally a compact landing page. Detailed operating procedures, mathematical assumptions, site-specific HPC instructions and validation evidence belong in this documentation tree.

When a solver change modifies equations, coefficients, boundary semantics, file formats or parallel ownership, update the corresponding documentation in the same pull request.
