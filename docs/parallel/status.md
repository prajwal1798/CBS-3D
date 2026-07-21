# Parallelisation status

## Current statement

CBS3D++ contains OpenMP acceleration, MPI initialisation, partition metadata, halo-exchange infrastructure and PETSc integration. The accepted production calculation is still serial/PETSc or shared-memory OpenMP.

The main MPI executable currently must not be described as a completed distributed CFD solver. Rank-local ownership, distributed assembly, reverse accumulation, forward halo updates and distributed pressure solution must all operate through CBS Steps 1–4 before that claim is valid.

## Implemented infrastructure

- optional MPI initialisation;
- MPI-aware build variants;
- partition metadata structures;
- owned/ghost-node concepts;
- halo-exchange implementation scaffolding;
- PETSc-enabled build path;
- DMPlex-oriented distributed-mesh development;
- OpenMP element colouring and threaded kernels;
- scaling-study and profiling utilities.

## Current blockers

1. The active main numerical path is not fully rank-local.
2. Mesh, field and output ownership are not yet consistently distributed through the full solver lifecycle.
3. The pressure system must use a distributed communicator and rank-owned rows/vectors.
4. Nodal assembly requires reverse additive communication from ghosts to owners.
5. Updated fields require forward owner-to-ghost communication between mathematical stages.
6. Boundary and material masks must be reconciled consistently across partition interfaces.
7. Parallel output must preserve global IDs and avoid duplicate ownership.
8. Serial and distributed results require controlled equivalence testing.

## Required communication pattern

For each field update, the distributed implementation should follow:

```text
forward halo of input field
-> assemble on owned elements
-> reverse additive halo of nodal contributions
-> update owned nodes and apply strong boundary conditions
-> forward halo of updated field
```

Communication belongs between mathematical stages. It should not change the CBS equation ordering.

## Ownership model

A rank-local partition requires explicit definitions for:

- owned elements;
- owned nodes/degrees of freedom;
- ghost nodes;
- local-to-global maps;
- neighbouring ranks;
- send/receive index lists;
- boundary-face ownership;
- material and phase ownership;
- global reductions for convergence and conservation.

Only owners may commit final nodal values. Ghost values are read-only replicas except during reverse additive assembly buffers.

## Development phases

### Phase 1 — distributed mesh acceptance

Acceptance criteria:

- global cell, vertex and boundary-face totals match the serial mesh;
- every global entity has consistent ownership;
- material and boundary labels survive distribution;
- ghost points and neighbour ranks are valid;
- one-, two- and four-rank diagnostics pass.

This phase does not execute CBS Steps 1–4.

### Phase 2 — distributed state and preprocessing

Acceptance criteria:

- rank-local state allocation;
- local geometry and shape derivatives;
- distributed boundary classification;
- consistent material-domain masks;
- local and global mass totals matching serial results.

### Phase 3 — distributed Step 1 and Step 3

Acceptance criteria:

- owned-element momentum assembly;
- reverse additive nodal accumulation;
- owner velocity updates;
- strong boundary enforcement;
- forward velocity halos;
- serial/distributed velocity equivalence.

### Phase 4 — distributed pressure solve

Acceptance criteria:

- distributed matrix/vector ownership;
- parallel PETSc communicator;
- correct pressure reference or null-space handling;
- rank-local RHS assembly;
- global KSP convergence reporting;
- pressure and divergence equivalence to the serial reference.

### Phase 5 — distributed Step 4 and turbulence

Acceptance criteria:

- fluid/solid thermal assembly across partitions;
- interface-node consistency;
- distributed turbulence-variable transport;
- global energy and scalar conservation.

### Phase 6 — parallel output and scaling

Acceptance criteria:

- correct parallel VTU/PVTU or equivalent output;
- deterministic global field reconstruction;
- strong-scaling results;
- weak-scaling results;
- load-balance and communication profiles;
- documented hardware, compiler, MPI, PETSc and partitioner.

## Verification requirements

Every distributed phase must compare against an accepted one-rank reference using:

```text
mesh and label totals
L2 and Linf field differences
mass-flow imbalance
pressure residual
velocity divergence
energy balance
iteration counts
convergence history
```

Floating-point operation order can produce small differences. Tolerances must be defined before the comparison; they must not be chosen after observing a failure.

## Scalability language

Permitted now:

> Designed for scalable distributed-memory execution; MPI domain decomposition is under active implementation.

Not permitted until Phase 6:

> Fully distributed solver.

> Demonstrated scalable parallel CFD solver.

> Linear or near-linear speed-up.

Published scalability claims must identify the case, mesh size, rank/thread layout, hardware, wall-clock region and baseline.
