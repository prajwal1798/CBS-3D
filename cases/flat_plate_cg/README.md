# Flat Plate - Native CG Comparison

Same 284,400-tetrahedron mesh and SA settings as `../flat_plate_petsc`, but
with `Pressure Field Solver = 1`.

This case exists only to measure the native CG baseline.  It stops at the
5,000-iteration cap without reaching the requested pressure tolerance.
