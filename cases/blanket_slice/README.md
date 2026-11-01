# Blanket slice: tapered plate, first wall, 4 nested U-channels

`BlanketSlice_CHT.geo` follows the design drawing (mm, x from the left back
corner):

| | Front-run top edge | Width | Left-leg outer edge at y = 0 | Fillet (both edges) |
|---|---|---|---|---|
| Plate | 280 | 350 (back), 301 (front), sides 85 deg | 0 | R6 |
| Channel 1 | 270 (first wall 10) | 15 | 20 | R10 |
| Channel 2 | 220 | 10 | 65 | R5 |
| Channel 3 | 150 | 10 | 105 | R5 |
| Channel 4 | 100 | 5 | 135 | R3.5 |

Channel legs are parallel to the plate sides; widths are true (perpendicular)
widths; the plate is symmetric about x = 175 mm. Through the thickness each
channel is square and centred in a 30 mm slice (7.5 mm walls either side of
channel 1). All values are `DefineConstant` parameters (Gmsh Parameters panel
or `-setnumber`).

IDs: volumes 10 fluid, 20 solid (21 armour if `tArm` > 0); surfaces 511
inlets (4 left legs), 520 outlets (4 right legs), 532 heated first wall (front
face + R6 fillets), 530 all other external faces.

Checked on the built CAD (`gmsh -0`, no mesh): every edge position, width and
fillet radius above is reproduced exactly; fluid 2.809e-4 m3, solid
2.453e-3 m3; inlet area 451.7 mm2 = (225 + 100 + 100 + 25) mm2 / sin 85 deg.

## Laminar study (stages L0, L1, L2)

| File | Content |
|---|---|
| `blanket_L0.par` | flow only (energy off) |
| `blanket_L1.par` | flow + energy, first-wall flux 1.0e5 W/m2 on BC 532 |
| `blanket_L2.par` | L1 + nuclear heating in the EUROFER, 1.6e6 exp(-d/0.15 m) W/m3, d = 0.28 m - y |
| `blanket_L1T.par` | L1 in physical time (global step, 60 s) for an unsteady flow; time averages via `CBS3D_AVERAGE_START_TIME` (runbook 6.4b) |
| `blanket.matprop` | 0 water at 298 K (997, 4182, 0.606, 8.9e-4), 1 EUROFER97 (7850, 500, 28.3) |

All stages: dimensional SI, mass-flow inlet BC 511 with 0.050 kg/s in total
at 298.15 K (the solver gives every inlet face the same normal speed, so
Re = 1870 / 1250 / 1250 / 620 in channels 1-4), outlets BC 520 at 0 Pa gauge,
steady stop on RelU,V,W < 1e-8 and RelT < 1e-9 after 20000 iterations. The
loads are one fifth of EU-DEMO nominal (0.5 MW/m2, ~8 MW/m3 at the first
wall): the largest that keep a laminar water wall below saturation in a
pressurised loop and EUROFER below 550 C. Expected: Q_in ~ 2.7 kW, coolant
rise ~ 12 K.

On Sunbird2 the whole chain is `jobs/sunbird/templates/blanket.sh`:

```bash
T=$CBS3D_CODE/jobs/sunbird/templates
bash $T/blanket.sh update ~/cbs3d-nuclear-heating.patch   # apply, build, self-check
bash $T/blanket.sh setup $CBS3D_SCRATCH/meshes/slice_h5.msh  # convert + partitions 40 80 128 256
bash $T/blanket.sh smoke L1                                # 300 iterations on 40 ranks
bash $T/blanket.sh run L1 --time 48:00:00                  # largest set that fits the free cores
bash $T/blanket.sh status
bash $T/blanket.sh check                                   # Q_out / Q_in, ports, residual plot
```

`blanket.sh` turns on the temperature limiter (`CBS3D_TEMPERATURE_AFC=1`,
algebraic flux correction, MPI runs only): with water (Pr ~ 6) the cell Peclet
number is several hundred, and on the check case below the unlimited scheme
undershot the 298.15 K inlet temperature by 21 K; with the limiter the minimum
stays at 298.15 K.

`CP_SCALE=5 bash $T/blanket.sh run L1` divides the EUROFER c_p by 5 for that
run only. The steady answer does not depend on the solid c_p; the solid then
reaches it about 5 times sooner in pseudo-time (the transient is not
physical). Converged = residuals flat AND `blanket.sh check` reports
Q_out / Q_in within 2 %.
