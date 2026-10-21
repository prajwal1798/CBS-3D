# Local Cases

The mesh bytes are not tracked.  The cases use hard links to the canonical
mesh in `C:\PhD\MESH\examples\3D\flat_plate` so there is exactly one copy on
disk.

## Recreate the links

```powershell
$source = "C:\PhD\MESH\examples\3D\flat_plate\flat_plate_sa_3d_corrected"
$cases = @{
    "flat_plate_local" = "flat_plate_local"
    "flat_plate_petsc" = "flat_plate_petsc"
}

foreach ($entry in $cases.GetEnumerator()) {
    $dir = Join-Path $PSScriptRoot $entry.Key
    $base = $entry.Value
    foreach ($ext in "plt", "bco", "material", "matprop") {
        $target = Join-Path $dir "$base.$ext"
        if (-not (Test-Path $target)) {
            New-Item -ItemType HardLink -Path $target -Target "$source.$ext" | Out-Null
        }
    }
}
```

If the mesh is on a different volume, replace `New-Item -ItemType HardLink`
with `Copy-Item`.

Validate a case:

```powershell
python tools\mesh_summary.py cases\flat_plate_local\flat_plate_local
```

## Cases

| Directory | Purpose |
| --- | --- |
| `flat_plate_local` | Native CG baseline, expected to stall |
| `flat_plate_petsc` | PETSc/AMG local production case |
| `flat_plate_p8` | 8-part mesh, MPI only |
