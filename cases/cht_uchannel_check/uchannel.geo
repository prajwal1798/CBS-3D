// ================================
// FILE    : uchannel.geo
// PURPOSE : SMALL CONJUGATE HEAT-TRANSFER CHECK CASE FOR THE BLANKET TOOLS:
//           A STEEL BLOCK (60 x 50 x 10 mm) COOLED BY ONE SQUARE U-CHANNEL
//           (5 mm) WHOSE INLET AND OUTLET BOTH OPEN ON THE FACE y = 0,
//           HEATED ON THE FACE y = H (BC 532) AND BY A VOLUMETRIC SOURCE
//           DECAYING FROM THAT FACE. AT STEADY STATE THE HEAT CARRIED OUT
//           BY THE WATER MUST EQUAL THE HEAT PUT IN (tools/energy_balance.py).
// ================================
SetFactory("OpenCASCADE");

DefineConstant[
    Wb = {0.060, Name "Check/block width (m)"},
    H  = {0.050, Name "Check/block height (m)"},
    Tz = {0.010, Name "Check/block thickness (m)"},
    a  = {0.005, Name "Check/channel side (m)"},
    xL = {0.015, Name "Check/left leg centre (m)"},
    xR = {0.045, Name "Check/right leg centre (m)"},
    yR = {0.035, Name "Check/channel run bottom (m)"},
    lc = {0.0015, Name "Check/element size (m)"}
];

z0 = 0.5 * (Tz - a);
Box(1) = {0, 0, 0, Wb, H, Tz};
Box(2) = {xL - a/2, 0, z0, a, yR + a, a};              // LEFT LEG
Box(3) = {xR - a/2, 0, z0, a, yR + a, a};              // RIGHT LEG
Box(4) = {xL - a/2, yR, z0, xR - xL + a, a, a};        // RUN
BooleanUnion(5) = { Volume{2}; Delete; }{ Volume{3, 4}; Delete; };
BooleanFragments{ Volume{1}; Delete; }{ Volume{5}; Delete; }

eps = 1e-6;
vFluid[] = Volume In BoundingBox{xL - a/2 - eps, -eps, z0 - eps, xR + a/2 + eps, yR + a + eps, z0 + a + eps};
vAll[] = Volume{:};
vSolid[] = vAll[];
vSolid[] -= vFluid[];

sIn[]  = Surface In BoundingBox{xL - a/2 - eps, -eps, z0 - eps, xL + a/2 + eps, eps, z0 + a + eps};
sOut[] = Surface In BoundingBox{xR - a/2 - eps, -eps, z0 - eps, xR + a/2 + eps, eps, z0 + a + eps};
sTop[] = Surface In BoundingBox{-eps, H - eps, -eps, Wb + eps, H + eps, Tz + eps};
bnd[] = Abs(CombinedBoundary{ Volume{vAll[]}; });
sWall[] = bnd[];
sWall[] -= sIn[];
sWall[] -= sOut[];
sWall[] -= sTop[];

Physical Volume("fluid", 10) = {vFluid[]};
Physical Volume("solid", 20) = {vSolid[]};
Physical Surface("inlet", 511) = {sIn[]};
Physical Surface("outlet", 520) = {sOut[]};
Physical Surface("heat_flux", 532) = {sTop[]};
Physical Surface("adiabatic_wall", 530) = {sWall[]};

Mesh.CharacteristicLengthMin = 0.5 * lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.MshFileVersion = 2.2;
