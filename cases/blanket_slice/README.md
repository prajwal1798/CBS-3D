# Blanket slice: tapered cooling plate with 4 nested U-channels

`BlanketSlice_CHT.geo`: a tapered plate (back 350 mm, front 270 mm, depth
400 mm, thickness 24 mm) with a first wall and armour on the front, cooled by
4 nested U-shaped square channels (12, 9.6, 7.8, 6 mm; largest next to the
first wall) with 12 mm gaps, ports on the back face, and the solid core inside
channel 4 removed (12 mm wall kept).

Uniform-offset construction: the plate outline, both edges of every channel
and the core cut-out are offsets of one rounded outline; every front-corner
arc shares one centre. Corners are concentric, channel widths and gaps are
constant through the bends. The plate fillet R0 must exceed the deepest
offset (checked; default R0 = 100 mm, deepest arc radius 10.6 mm).

IDs: volumes 10 fluid, 20 structure, 21 armour; surfaces 511 inlets (4),
520 outlets (4), 532 heated front (front face + fillets, 5 faces), 530 other.

Checked (geometry only, `gmsh -0`): 4 fluid volumes (3.005e-4 m3), structure
1.550e-3 m3, armour 6.80e-6 m3; inlet area 3.347e-4 m2 = sum of the four
squares x 1/cos(taper).
