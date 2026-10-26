# Box - pipeline smoke case

Unit cube, 1,089 tetrahedra. Proves the whole chain in seconds:
Gmsh CLI -> CBS converter -> mesh check -> serial, PETSc and MPI runs.

```bash
cd ~/PhD/CFD/CBS-3D
bash tools/cbs_case.sh cases/box/box.geo cases/box/box.par cases/box/box.matprop \
    ~/PhD/runs/box --parts 4 --python ~/venvs/phd/bin/python
```

Expected: `1089 tets, 291 nodes, 400 boundary faces`,
`BC face counts: 511=68, 520=66, 530=266`, `RESULT: PASS`.
Gmsh 4.12-4.15 give 1089; some versions give 1091, which is fine.

`box.par` uses the native CG solver (`Pressure Field Solver` = `1`).
Change that line to `3` for PETSc. After 5 iterations every build gives
`RelU=5.214855e-03  RelV=1.202367e-01  RelW=1.732350e-01  RelP=4.308202e-01`.

See `docs/LINUX_SUNBIRD_RUNBOOK.md`.
