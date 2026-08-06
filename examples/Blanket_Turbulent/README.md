# Blanket turbulent CHT example

## Status

`DEVELOPMENT ONLY`

This case is reserved for the Spalart-Allmaras blanket calculation. It must not
be advertised as a supported distributed example until the MPI SA transport,
halo exchange, owner constraints, convergence metrics, restart data and output
fields have been implemented and validated.

## Current execution restriction

- serial SA development: permitted for controlled validation;
- MPI SA production: explicitly rejected by the parallel application;
- hybrid MPI/OpenMP SA: unsupported.

## Required validation stages

1. zero-eddy-viscosity laminar recovery;
2. wall-distance verification;
3. manufactured or canonical SA transport check;
4. flat-plate skin-friction comparison;
5. serial mesh-refinement study;
6. serial versus MPI rank-count agreement;
7. turbulent CHT energy-balance validation;
8. long-run restart/checkpoint validation.

Authoritative inputs will be promoted only after the scratch inventory is
complete.
