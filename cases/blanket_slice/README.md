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
