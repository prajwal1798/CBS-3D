# Developer architecture

## Source layout

```text
include/cbs/core/         Fundamental types, state and run configuration
include/cbs/io/           Mesh input and post-processing interfaces
include/cbs/preprocess/   Geometry and finite-element preprocessing
include/cbs/boundary/     Boundary classification and enforcement
include/cbs/timestep/     Time-step and diagonal updates
include/cbs/assembly/     Momentum, pressure and energy assembly
include/cbs/linalg/       Native and PETSc linear algebra
include/cbs/solver/       CBS steps, convergence and solver orchestration
include/cbs/parallel/     Partition metadata and halo exchange
include/cbs/utils/        Profiling and supporting utilities

src/                      Implementations following the same module structure
tools/                    Mesh, monitoring, post-processing and scaling tools
docs/                     User, numerical, validation and development records
```

## Active call order

```text
main
└── Solver::run
    ├── Solver::initialise
    │   ├── MeshIO::readAll
    │   ├── boundary and input validation
    │   ├── shape-function derivatives and geometry
    │   ├── boundary-face assignment and normals
    │   ├── lumped mass and element-size preprocessing
    │   ├── material and domain classification
    │   ├── element colouring
    │   └── initial boundary enforcement
    ├── Solver::preparePressureSystem
    └── iteration loop
        ├── copy iteration-old fields
        ├── compute/update time-step diagonals
        ├── Steps::step1
        ├── Steps::step2
        ├── Steps::step3
        ├── Steps::step4, when enabled
        ├── update derived fields
        ├── evaluate convergence
        └── write residuals and solution output
```

The mathematical order is part of the formulation. Refactoring must not move communication, boundary enforcement or field overwrites across stage boundaries without a numerical review.

## State model

The principal state contains:

- mesh connectivity and coordinates;
- boundary records and mappings;
- velocity, pressure and temperature fields;
- iteration-old fields;
- intended physical-time history arrays;
- momentum and scalar RHS storage;
- tetrahedral gradients and determinants;
- lumped mass and time diagonals;
- element/global pressure coefficients;
- material and domain masks;
- element-colouring groups;
- MPI context and developing partition metadata.

The implementation deliberately preserves one-based indexing in major scientific arrays. Index zero is unused. A conversion to zero-based indexing is a numerical migration requiring complete regression, not a style-only change.

## Assembly discipline

### OpenMP

Shared-node element scatter is protected by element colouring. Kernels operating on colour groups may update nodal arrays without atomics when no two elements in a colour share a node.

Any new assembly loop must document:

```text
read fields
written arrays
scatter pattern
race-prevention method
material-domain restriction
boundary contribution
```

### MPI

Distributed assembly must separate:

1. local owned-element computation;
2. reverse additive accumulation to owners;
3. owner-only nodal update;
4. strong boundary enforcement;
5. forward owner-to-ghost update.

See [parallel status](../parallel/status.md).

## Linear algebra

### Native pressure path

- compact element pressure coefficients;
- matrix-free operator application;
- preconditioned conjugate gradient;
- global pressure boundary treatment.

### PETSc pressure path

- PETSc system setup and KSP solution;
- runtime preconditioner selection;
- HYPRE support when available.

The current PETSc implementation must not cache an operator across coefficient changes unless the cache key includes all quantities affecting the matrix, including time-step scaling.

## Boundary contract

Boundary behaviour spans input, preprocessing and solver steps. A boundary identifier is not fully implemented unless it is:

1. accepted by the input reader;
2. mapped by `.bco`;
3. classified during preprocessing;
4. assembled where required;
5. enforced at the correct CBS stage;
6. represented correctly in output/debug fields;
7. covered by a regression case.

## Material contract

The active convention is:

```text
material ID 0      fluid
non-zero ID        solid
```

The convention is embedded in domain masks and must be generalised before multiple independently identified fluid materials are supported.

## Adding a numerical feature

A numerical pull request should contain:

1. governing equation or algorithm reference;
2. discrete form and coefficient definitions;
3. affected files and arrays;
4. units/dimensional assumptions;
5. boundary and material behaviour;
6. serial regression;
7. OpenMP race analysis;
8. MPI communication requirements;
9. validation or verification case;
10. documentation update.

## Dormant and unsupported paths

Configuration parsing does not guarantee implementation. Controls that are ignored, rejected later or only partially connected must be marked unsupported until the full call path and validation exist.

Do not advertise dormant banded pressure, restart, true transient, explicit CBS or partially connected stabilisation options merely because storage or parser fields are present.

## Change review checklist

- Does the code compile with warnings enabled?
- Did a state array change ownership or lifetime?
- Did the mathematical stage order change?
- Are all nodal scatters race-free?
- Are fluid/solid masks correct?
- Are units consistent?
- Are input files still backward compatible?
- Does serial regression pass?
- Is the distributed implication documented?
- Are README/status claims still accurate?
