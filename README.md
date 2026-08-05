## Problem

BC511 mass-flow inlet velocity was normalised using the complete geometric
inlet area before final nodal boundary-condition priority was resolved.

Some inlet rim nodes were subsequently overridden by solid/interface/wall
conditions and assigned zero velocity. Therefore, the actual P1-integrated
mass flow was lower than the requested value.

For the straight-pipe CHT benchmark:

- Requested mass flow: 1.000000e-03 kg/s
- Previous delivered mass flow: approximately 9.18e-04 kg/s
- Error: approximately 8.2%

## Correction

This change:

- resolves final velocity boundary-condition priority before mass-flow scaling;
- identifies BC511 nodes that remain active after priority resolution;
- assembles the actual discrete P1 inlet support;
- includes any fixed higher-priority velocity contribution;
- computes the required profile amplitude from the discrete flux equation;
- stores and applies the final nodal BC511 velocity consistently;
- independently integrates the resulting P1 flux for verification;
- supports both serial and distributed preprocessing paths.

No mesh-dependent empirical correction factor is used.

## Validation

Windows Release build:

- Visual Studio 2022
- C++20
- OpenMP enabled
- MPI disabled
- PETSc disabled

One-iteration validation:

- BC511 inlet faces: 655
- Geometric inlet area: 3.136387096911e-04
- Active P1 support area: 2.880603406112e-04
- Topological inlet nodes: 360
- Active inlet nodes: 297
- Overridden inlet nodes: 63
- Computed profile amplitude: 3.481940710593e-03
- Requested mass flow: 1.000000000000e-03 kg/s
- Integrated P1 mass flow: 1.000000000000e-03 kg/s
- Relative verification error: 3.686287386451e-15



## Remaining validation

The MPI/PETSc distributed path must be compiled and tested on Sunbird before
this PR is merged.
