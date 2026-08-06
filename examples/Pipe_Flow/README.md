# Straight-pipe laminar CHT benchmark

## Status

`MPI/PETSc PRODUCTION VALIDATED`

This example is the current authoritative validation case for the corrected
BC511 discrete mass-flow normalisation and the distributed laminar CHT path.

## Validated production configuration

- global mesh: 648,353 tetrahedra and 119,823 nodes;
- MPI ranks: 40;
- compute nodes: 5;
- MPI ranks per node: 8;
- OpenMP threads per rank: 1;
- fixed pseudo-time step: `1.0e-3`;
- iterations: 1,000,000;
- measured elapsed time: 5 h 06 min 28 s;
- requested inlet mass flow: `1.0e-3 kg/s`.

## BC511 discrete verification

- inlet faces: 655;
- geometric inlet area: `3.136387096911e-04`;
- active P1 support area: `2.880603406112e-04`;
- topological inlet nodes: 360;
- active profile nodes: 297;
- overridden inlet nodes: 63;
- integrated P1 mass flow: `1.000000000000e-03 kg/s`;
- relative verification error: `3.686287386451e-15`.

## Hydrodynamic reference checks

- axial pressure follows the Hagen-Poiseuille gradient;
- numerical inlet gauge pressure: approximately `0.26458 Pa`;
- analytical inlet gauge pressure: approximately `0.25541 Pa`;
- inlet pressure difference: approximately `3.59%`;
- radial velocity profile closely follows
  `u/Ubar = 2[1 - (r/R)^2]`.

## Pressure units

The solver stores incompressible kinematic pressure internally. Dimensional
post-processing must use:

```text
p_gauge_Pa = density * p_kinematic
```

For this case, the inlet density is `997 kg/m^3`.

## Remaining packaging work

The authoritative serial input set, partition manifest, Slurm template,
post-processing scripts, accepted CSV values and file checksums will be copied
from the retained corrected production run after the Sunbird workspace
inventory identifies the final source directory.
