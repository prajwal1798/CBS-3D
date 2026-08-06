# Blanket laminar CHT example

## Status

`IMPORT PENDING`

This directory will contain the authoritative laminar blanket input set after
the Sunbird workspace inventory identifies the correct mesh, material files,
boundary-condition files, partition set and validated output manifest.

## Intended execution modes

- serial reference;
- serial OpenMP;
- serial PETSc pressure solve;
- MPI/PETSc distributed production;
- hybrid MPI/OpenMP after scaling validation.

## Promotion requirements

Before this case is marked production-ready, record:

- source geometry and mesh-generation commit;
- global tetrahedron, node and boundary-face counts;
- fluid/solid material counts;
- all boundary-condition identifiers;
- inlet mass flow, outlet pressure and thermal loading;
- serial/MPI agreement;
- mass and energy balance;
- restart/checkpoint validation;
- checksums for all large external assets.

Large partitions and result fields should be referenced from `manifest.json`
rather than duplicated in Git history.
