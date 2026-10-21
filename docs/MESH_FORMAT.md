# Native CBS3D Mesh Format

A case is a set of files that share one base name:

```
<case>.plt        tetrahedra, node coordinates, boundary faces
<case>.bco        declared boundary-condition codes
<case>.material   per-element material id (required when material mode is on)
<case>.matprop    material properties
<case>.par        solver controls
<case>.mpi        optional rank-local partition metadata
```

All indices are **one-based**.

## `<case>.plt`

The file is plain text with three sections.

### Header

```
nelem npoin nboun
```

### Tetrahedra

```
id n1 n2 n3 n4
```

`id` runs from 1 to `nelem`.  The four node ids are the vertices of a linear
tetrahedron.

### Nodes

```
id x y z
```

The production writer uses fixed-width fields.  A negative x coordinate can
touch the node id, for example:

```
         1-5.0000000000000000e-01 0.0000000000000000e+00 0.0000000000000000e+00
```

Parsers must therefore not assume a space between the id and the first
coordinate.  `tools/mesh_summary.py` handles this form explicitly.

### Boundary faces

```
n1 n2 n3 parent_element bc_code
```

The face is a triangular face of `parent_element`.  Boundary faces are listed
once.  Interior faces are implied by the tetrahedra and are not listed.

## `<case>.bco`

The first line is:

```
nflag mass_flow_flag
```

The next `nflag` lines declare the active boundary-condition codes:

```
code code
```

Example:

```
         4         0
       506       506
       511       511
       520       520
       530       530
```

The `.bco` file declares which codes are active.  The actual faces and their
parent elements are in the `.plt` file.

## Boundary-condition codes used by the flat plate

| Code | Meaning |
| --- | --- |
| `500` | Adiabatic wall with prescribed velocity |
| `501` | Constant-temperature wall |
| `502` | Constant-temperature wall (second variant) |
| `506` | Symmetry |
| `510` | Prescribed inlet velocity and temperature |
| `511` | Mass-flow inlet |
| `520` | Pressure outlet |
| `530` | No-slip adiabatic wall |
| `901` | Conformal fluid/solid CHT interface |

The flat plate uses `506` symmetry, `511` mass-flow inlet, `520` pressure
outlet and `530` no-slip adiabatic wall.

The canonical flat-plate controls are `Re = 1.0e6`, `U = 1`, `T = 300`,
`rho = 1` and `mu = 1.0e-6`.  The free-stream line is:

```
Free stream values - Ux Uy Uz P T
1.0  0.0  0.0  0.0  300.0
```

Keep this consistent with the prescribed inlet velocity `1.0 0.0 0.0`.  A
free-stream value of `20` from an older case makes the physical Reynolds
number twenty times the intended value and changes the SA near-wall behaviour.

## `<case>.material`

One record per tetrahedron:

```
element_id n1 n2 n3 n4 material_id
```

`material_id` is zero for fluid.  Solid ids index the `.matprop` table.

## `<case>.matprop`

The first line is the number of materials.  Each following line is:

```
id name phase density cp conductivity ... 
```

Example:

```
2
0 Air_Unit fluid 1.0 1.0 1.0 1e-06 0.0
1 RAFM_EUROFER97 solid 7850.0 500.0 28.3 0.0 0.0
```

## `<case>.mpi`

The `.mpi` file is present for rank-local partitions.  It is not a substitute
for a real MPI run: a rank-local mesh contains only that rank's elements and
owned plus ghost nodes.

Example header:

```
CBS3D_MPI_PARTITION 1
PARTITION 1 0 8
GLOBAL 284400 64320 66548
LOCAL 35363 8268 8268 0 8516
NODES
1 6 0
2 12 0
...
```

- `PARTITION <id> <rank> <parts>`
- `GLOBAL <tetrahedra> <nodes> <boundary_faces>`
- `LOCAL <tetrahedra> <nodes> <owned_nodes> <ghost_nodes> <boundary_faces>`
- `NODES` starts the node-owner/ghost table.

## Validation

`tools/mesh_summary.py` performs a Gmsh-free integrity check:

```powershell
python tools\mesh_summary.py cases\flat_plate_local\flat_plate_local `
    --json cases\flat_plate_local\mesh_summary.json
```

It checks:

- header counts against the actual number of records;
- one-based, contiguous node ids;
- tetrahedron connectivity and positive volume;
- face manifoldness (no face shared by more than two elements);
- every boundary face is a face of its declared parent element;
- `.bco` codes match the boundary faces in `.plt`;
- `.material` element ids and material ids;
- `.mpi` local counts against the `.plt` counts.

The 284,400-cell flat plate reports:

```
Tetrahedra    : 284400
Nodes         : 64320
Boundary faces: 66548
Total volume  : 5.000000000000e-02
Bounding box  : x[-0.5, 2] y[0, 1] z[0, 0.02]
BC face counts: 506=64640, 511=474, 520=474, 530=960
Material tets : 0=284400
RESULT: PASS
```

## Partition workflow

The canonical mesh examples are in `C:\PhD\MESH\examples\3D\flat_plate`.
Gmsh itself is not installed as a CLI on the tower; the Python `gmsh` module
(4.15.2) is used instead.

### 1. Partition the mesh

```powershell
python C:\PhD\MESH\examples\3D\flat_plate\tools\partition_mesh.py `
    C:\PhD\MESH\examples\3D\flat_plate\flat_plate_sa_3d_corrected.msh `
    --parts 8 `
    --out-dir C:\PhD\CFD\CBS3D_MASTER\cases\flat_plate_p8\partitioned_8 `
    --name flat_plate_p8
```

The tool uses METIS through the Gmsh Python module.  The 8-part flat plate
has about 7,157 edge cuts and a load imbalance of about 1.018.

### 2. Export rank-local CBS cases

```powershell
python C:\PhD\MESH\examples\3D\flat_plate\tools\export_rank_local_cbs.py `
    C:\PhD\CFD\CBS3D_MASTER\cases\flat_plate_p8\partitioned_8\flat_plate_p8_partitioned_8.msh `
    --parts 8 `
    --case-name flat_plate_p8 `
    --out-dir C:\PhD\CFD\CBS3D_MASTER\cases\flat_plate_p8\ranks `
    --par C:\PhD\MESH\examples\3D\flat_plate\flat_plate_sa_3d_corrected.par `
    --matprop C:\PhD\CFD\CBS3D_MASTER\cases\flat_plate_local\flat_plate_local.matprop
```

Validate each rank:

```powershell
python tools\mesh_summary.py `
    cases\flat_plate_p8\ranks\rank_0000\flat_plate_p8_rank_0000
```

Rank 0 reports 35,363 tetrahedra, 8,268 nodes and 8,516 boundary faces with
`RESULT: PASS`.

> A rank-local case must be run with MPI.  Running one rank-local file as a
> serial case silently omits all other partitions and their halo exchanges.
