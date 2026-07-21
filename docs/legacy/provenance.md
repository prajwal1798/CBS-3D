# Legacy provenance

## Summary

CBS3D++ is a modern C++ redevelopment in the numerical lineage of the legacy CBSFlow 3D Fortran solver. The project retains the Characteristic-Based Split finite-element methodology while replacing the software architecture, build system, data structures, linear-algebra interfaces, post-processing and parallel-development infrastructure.

This page records lineage. It does not transfer or expand rights in any legacy material.

## Public legacy description

Zeta Computational Resources publicly describes CBSFlow 3D as:

- a finite-element solver for the Navier–Stokes equations;
- capable of incompressible-flow simulation;
- based on characteristic-based split convection stabilisation;
- using a fractional-step approach;
- available historically in serial and MPI-oriented forms;
- associated with linear tetrahedral-element benchmark data.

Primary public references:

- [CBSFlow 3D parallel description](https://www.zetacomp.com/software/cbsflow3d.asp)
- [CBSFlow 3D serial description](https://www.zetacomp.com/software/cbsflow3dserial.asp)
- [Zeta software page and terms](https://www.zetacomp.com/software.asp)
- [Zeta 3D benchmark data](https://www.zetacomp.com/benchmarks/benchmarks3d.asp)

Methodological reference:

O. C. Zienkiewicz, R. L. Taylor and P. Nithiarasu, *The Finite Element Method for Fluid Dynamics*, 7th edition.

## Copyright boundary

Zeta states that its source codes and associated ideas are copyrighted, that redistribution of source is not permitted, and that reproduction requires written permission.

Accordingly:

- legacy source must not be added to this public repository unless redistribution permission is documented;
- downloaded legacy manuals and case archives must not be republished automatically;
- public links and bibliographic references should be used where possible;
- newly generated regression summaries should not expose protected source text;
- the repository licence must not imply relicensing of third-party legacy material.

## Required rights review

Before the first formal public release, the project owner should record:

1. how the legacy source was obtained;
2. the licence or written permission attached to it;
3. whether the C++ implementation is a translation, derivative work or independently written implementation;
4. which files, if any, contain translated legacy expression;
5. whether those files may be publicly distributed;
6. the copyright holders and required notices;
7. the licence applicable to newly authored code.

Until that review is complete, no broad open-source licence should be attached without legal and institutional confirmation.

## Numerical continuity

The redevelopment intends to preserve the accepted behaviour of the legacy semi-implicit incompressible method, including:

- four-stage CBS update order;
- finite-element shape-function conventions;
- pressure-correction semantics;
- boundary-condition interpretation;
- legacy case-file compatibility where practical;
- benchmark comparability.

Continuity must be demonstrated through the [legacy regression protocol](../validation/legacy-regression.md), not inferred from similar function names.

## Modern implementation scope

The modern project introduces or restructures:

- C++20 source organisation;
- typed state containers;
- modular mesh, boundary, assembly and solver components;
- CMake and Make build variants;
- OpenMP element assembly;
- PETSc KSP/HYPRE integration;
- MPI ownership and halo-exchange infrastructure;
- multiphysics material handling;
- turbulence-model development;
- VTU/PVD output;
- residual monitoring and profiling;
- controlled validation and HPC procedures.

## Intentional-difference register

The following items require explicit tracking because they can differ from a legacy implementation:

| Area | Required record |
|---|---|
| Momentum stress form | Exact active viscous operator |
| Higher-order CBS terms | Terms retained or omitted in each step |
| Time integration | Pseudo-time versus true transient formulation |
| Pressure solver | Native CG versus PETSc operator and tolerance |
| Material classification | Fluid/solid material-ID convention |
| Boundary semantics | Mapping and enforcement order |
| Parallel reductions | Ownership and accumulation order |
| Output | Field naming, indexing and sampling |

## Attribution statement

A release-ready attribution statement should distinguish clearly between:

- the CBS method and its published authors;
- the legacy CBSFlow software and its copyright holder;
- the CBS3D++ redevelopment and its contributors;
- third-party libraries such as PETSc, HYPRE, MPI, Gmsh and ParaView.

This statement will be finalised with the repository licence and institutional approval.
