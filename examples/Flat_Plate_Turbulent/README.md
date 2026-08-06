# Turbulent flat-plate benchmark

## Status

`SERIAL SPALART-ALLMARAS DEVELOPMENT`

This example is the primary canonical case for validating the current
Spalart-Allmaras implementation before turbulence is admitted to the MPI
production loop.

## Current execution modes

- serial single-thread: development validation;
- serial OpenMP: development validation;
- serial PETSc pressure solve: to be regression-tested;
- MPI and hybrid SA: unsupported and explicitly rejected.

## Required reference quantities

The promoted case must report at minimum:

- freestream velocity, density and viscosity;
- plate length and Reynolds-number range;
- wall-distance accuracy;
- `nu_tilde`, turbulent viscosity and SA source-term bounds;
- local and integrated skin-friction coefficient;
- comparison with the selected laminar/turbulent flat-plate correlation;
- velocity profiles in wall coordinates where appropriate;
- mesh-spacing and near-wall resolution study;
- residual and global conservation histories.

The exact mesh and input set will be selected from the retained development
folders after the non-destructive Sunbird inventory is complete.
