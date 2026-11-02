// G4: conformal two-material heat-conduction benchmark.
// Fluid at 0 <= y <= Lf; solid at Lf <= y <= Lf+Ls.
// Zero imposed velocity. Temperature T=0 at y=0 (BC 502), inward
// Neumann heat flux q at y=Lf+Ls (BC 532), remaining external
// surfaces zero normal conductive flux (BC 530 except pressure BC 520).
// The internal y=Lf interface is NOT a physical boundary:
// conformal tetrahedra share nodes and the weak form enforces flux continuity.
// Exact stationary field (for zero advection and zero volumetric heating):
// T(y) = q*y/kf                         for 0 <= y <= Lf
//      = q*Lf/kf + q*(y-Lf)/ks          for Lf < y <= Lf+Ls.
SetFactory("OpenCASCADE");
DefineConstant[
 W  = {0.004, Name "Slab/width"},
 Lf = {0.001, Name "Slab/fluid thickness"},
 Ls = {0.003, Name "Slab/solid thickness"},
 lc = {0.00075, Name "Slab/mesh size"}
];
Box(1) = {0, 0, 0, W, Lf, W};
Box(2) = {0, Lf, 0, W, Ls, W};
// Fragmentation makes the fluid/solid interface conformal.
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }

eps = 1e-7;
fluid[] = Volume In BoundingBox{-eps, -eps, -eps, W+eps, Lf+eps, W+eps};
solid[] = Volume In BoundingBox{-eps, Lf-eps, -eps, W+eps, Lf+Ls+eps, W+eps};

cold[] = Surface In BoundingBox{-eps, -eps, -eps, W+eps, eps, W+eps};
hot[]  = Surface In BoundingBox{-eps, Lf+Ls-eps, -eps, W+eps, Lf+Ls+eps, W+eps};
iface[] = Surface In BoundingBox{-eps, Lf-eps, -eps, W+eps, Lf+eps, W+eps};
press[] = Surface In BoundingBox{W-eps, -eps, -eps, W+eps, Lf+eps, W+eps};

all[] = Surface{:};
walls[] = all[];
walls[] -= cold[];
walls[] -= hot[];
walls[] -= iface[];
walls[] -= press[];

Physical Volume("fluid", 10) = {fluid[]};
Physical Volume("solid", 20) = {solid[]};
Physical Surface("cold_wall", 502) = {cold[]};
Physical Surface("heat_flux", 532) = {hot[]};
Physical Surface("outlet", 520) = {press[]};
Physical Surface("adiabatic_wall", 530) = {walls[]};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.MshFileVersion = 4.1;
