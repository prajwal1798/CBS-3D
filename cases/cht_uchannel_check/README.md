# Conjugate U-channel check case

A small stand-in for the blanket slice that runs in minutes: a EUROFER block
60 x 50 x 10 mm cooled by one square 5 mm U-channel whose inlet and outlet
both open on y = 0, heated on y = H by 1e5 W/m2 (BC 532) and by
1.6e6 exp(-(H - y)/0.02 m) W/m3 in the steel (the blanket L2 `.par` with the
geometry values changed). Water 1.25e-3 kg/s at 298.15 K (Re ~ 250).

It exercises everything the blanket run uses: mass-flow inlet, pressure
outlet, conformal fluid-solid interface, wall flux, the nuclear-heating block,
the temperature limiter, MPI partitions, and `tools/energy_balance.py`.

```bash
gmsh uchannel.geo -3 -format msh2 -o uchannel.msh
python3 ../../tools/gmsh_to_cbs.py uchannel.msh --case-name uchannel --out-dir serial \
    --par uchannel.par --matprop uchannel.matprop
python3 ../../tools/export_rank_local_cbs.py serial/uchannel --parts 2 --out-dir p2 --gmsh gmsh
mkdir run && cd run
CBS3D_TEMPERATURE_AFC=1 mpirun -np 2 cbs3d_parallel uchannel ../p2
python3 ../../../tools/energy_balance.py output/uchannel --log solver.out
```

## Measured (45 381 tets, 2 MPI ranks, solid c_p / 10, temperature AFC on)

| Iterations | Q_out / Q_in | T range [K] |
|---|---|---|
| 10 000 | 0.453 | 298.15 - 375.1 |
| 20 000 | 0.774 | 298.15 - 395.0 |
| 30 000 | 0.911 | 298.15 - 402.5 |
| 40 000 | 0.971 | 298.15 - 405.5 |
| 50 000 | 0.996 | 298.15 - 406.8 |

Q_in = 76.25 W (60 W flux + 16.25 W nuclear, printed by the solver). The
heat carried out by the water converges to the heat put in, so flux, source,
interface coupling and limiter are conservative.

Two findings that matter for the blanket:

- Without the limiter (serial run, same case) the water temperature fell to
  277.2 K, 21 K below the 298.15 K inlet: a convection-dominated undershoot
  (cell Peclet ~ 500). With `CBS3D_TEMPERATURE_AFC=1` the minimum stays
  298.15 K. `blanket.sh` turns the limiter on.
- The outlet mass flow is 2.7 % below the inlet (constant from 10 000
  iterations on): the discrete continuity error of a mesh with only 3-4
  cells across the 5 mm channel. `energy_balance.py` reports it separately;
  on the blanket meshes (0.75 mm wall cells) it should be smaller, and it is
  a mesh-refinement indicator.
