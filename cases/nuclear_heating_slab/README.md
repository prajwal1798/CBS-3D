# Nuclear-heating verification slab

A EUROFER bar, 0 <= y <= L (L = 0.10 m, 10 x 10 mm section), with

- heat flux q = 1e5 W/m2 entering at y = L (BC 532, the "first wall"),
- T = 0 at y = 0 (BC 502), all other faces adiabatic (BC 530),
- volumetric source Q = A exp(-(L - y)/lambda), A = 1.6e6 W/m3, lambda = 0.05 m
  (optional `.par` block, depth measured from y = L along -y).

A separate stagnant water cube only gives the solver its pressure reference.
The exact steady temperature is

    T(y) = (A/k) [lambda y - lambda^2 (exp(-(L-y)/lambda) - exp(-L/lambda))] + q y / k

with T(L) = 513.83 K for k = 28.3 W/m K, and the heat input is
A W^2 lambda (1 - exp(-L/lambda)) + q W^2 = 16.918 W.

Measured (lc = 5 mm, 544 tets, steady stop at 89 434 iterations):

| Run | max abs(T - T_exact) | printed heat input |
|---|---|---|
| serial | 0.21 K (4.2e-4 of T(L)) | 16.916 W (1e-4) |
| 2 MPI ranks | same field to 2.3e-13 K | 16.916 W |
| 2 MPI ranks, `CBS3D_TEMPERATURE_AFC=1` | same field to 2.3e-13 K (AFC limits fluid elements only) | 16.916 W |
| 2 MPI ranks, AFC on solids too (`CBS3D_AFC_SOLID=1`, the behaviour before 2026-11-01) | 3.3 K off serial, stops at 88 597 | 16.916 W |

`jobs/sunbird/templates/selfcheck_slab.slurm` repeats this (serial and 4
ranks) and prints `SELF-CHECK PASS`; `tools/check_nuclear_slab.py` does the
comparison.
