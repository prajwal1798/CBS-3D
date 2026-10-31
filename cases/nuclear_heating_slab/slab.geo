// ================================
// FILE    : slab.geo
// PURPOSE : VERIFICATION OF THE DEPTH-DEPENDENT VOLUMETRIC SOURCE.
//           ONE SOLID (EUROFER) BAR, 0 <= y <= L.  HEAT FLUX q ENTERS AT
//           y = L (THE "FIRST WALL", BC 532), T = 0 AT y = 0 (BC 502),
//           ALL OTHER FACES ADIABATIC (BC 530).  SOURCE Q = A exp(-(L-y)/lam).
//           EXACT STEADY SOLUTION (1D):
//             T(y) = (A/k)[lam y - lam^2 (exp(-(L-y)/lam) - exp(-L/lam))] + q y / k
//           A SEPARATE STAGNANT FLUID CUBE (NOT TOUCHING THE BAR) ONLY GIVES THE
//           SOLVER THE PRESSURE REFERENCE IT NEEDS (BC 520 ON ITS TOP FACE).
// ================================
SetFactory("OpenCASCADE");

DefineConstant[
    L   = {0.10,  Name "Slab/L (m)"},
    W   = {0.01,  Name "Slab/cross-section side (m)"},
    lc  = {0.005, Name "Slab/element size (m)"}
];

Box(1) = {0, 0, 0, W, L, W};
Box(2) = {2*W, 0, 0, W, W, W};

eps = 1e-5;
yLow[]  = Surface In BoundingBox{-eps, -eps, -eps, W+eps, eps, W+eps};
yHigh[] = Surface In BoundingBox{-eps, L-eps, -eps, W+eps, L+eps, W+eps};
fTop[]  = Surface In BoundingBox{2*W-eps, W-eps, -eps, 3*W+eps, W+eps, W+eps};
all[]   = Surface{:};
walls[] = all[];
walls[] -= yLow[];
walls[] -= yHigh[];
walls[] -= fTop[];

Physical Volume("solid", 20)      = {1};
Physical Volume("fluid", 10)      = {2};
Physical Surface("cold_wall", 502) = {yLow[]};
Physical Surface("heat_flux", 532) = {yHigh[]};
Physical Surface("outlet", 520)    = {fTop[]};
Physical Surface("adiabatic_wall", 530) = {walls[]};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.MshFileVersion = 4.1;
