# Flat Plate - 8-Part MPI Mesh

Eight METIS partitions of the full flat plate.

- `partitioned_8/` contains the partitioned Gmsh mesh.
- `ranks/rank_XXXX/` contains the exported rank-local CBS files.
- Global mesh: 284,400 tetrahedra, 64,320 nodes, 66,548 boundary faces.
- Rank 0: 35,363 tetrahedra, 8,268 nodes, 8,516 boundary faces.
- Edge cuts: about 7,157; load imbalance about 1.018.

Validation:

```powershell
python tools\mesh_summary.py cases\flat_plate_p8\ranks\rank_0000\flat_plate_p8_rank_0000
```

This case must be run with MPI.  A single rank-local file is not a valid
serial case; it omits the other partitions and their halo exchanges.  The
local MSYS2 PETSc is an MPIUNI build, so use Sunbird or a local MPI PETSc
build for this case.

The mesh bytes and rank directories are generated artefacts and are not
tracked.  Regeneration commands are in `docs/MESH_FORMAT.md`.
