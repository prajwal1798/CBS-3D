//==============================================================================
// ManifoldPanel_CHT.geo
//
// PARAMETRIC MANIFOLD-FED COOLED PANEL: A CONJUGATE-HEAT-TRANSFER BENCHMARK
// GEOMETRY FOR CBS3D++_SI.
//
// A layered solid panel, heated on its front face, cooled by N parallel square
// channels that are fed by an inlet header and drained by an outlet header.
// The same arrangement appears in actively cooled walls (fusion first-wall and
// blanket panels), electronics and battery cold plates and heat-exchanger
// plates. It combines four effects a single-channel cell cannot show:
//   1. FLOW DISTRIBUTION between parallel channels (header pressure losses);
//   2. LATERAL HEAT SPREADING between neighbouring channels through the solid;
//   3. LAYERED SOLIDS (optional armour layer on the heated face);
//   4. U-TYPE OR Z-TYPE MANIFOLDS (ports on the same or on opposite sides).
//
// PLAN VIEW (z = channel mid-height), U-TYPE (manifold = 0):
//
//      y ^      inlet header        N channels (pitch p)       outlet header
//        |     +------+-------------------------------------+------+
//        |     |      |=====================================|      |
//        |     |  ->  |=====================================|  ->  |
//        |     |      |=====================================|      |
//        |     |      |=====================================|      |
//        |     +--  --+-------------------------------------+--  --+
//        +---> x  inlet port (y = -Wp/2)          outlet port: y = -Wp/2 (U)
//                                                              y = +Wp/2 (Z)
// SECTION (x = const):
//        +------------------------------------------------+  z = Ttop
//        |   armour  (tArm, optional)  heated face above   |
//        +------------------------------------------------+
//        |   front wall (tFront)                           |
//        |   [a]   [a]   [a]   [a]   ...  channels         |
//        |   back wall (tBack)                             |
//        +------------------------------------------------+  z = 0
//
// ALL LENGTHS ARE MULTIPLES OF THE CHANNEL SIDE a (= HYDRAULIC DIAMETER), SO
// THE BENCHMARK IS DEFINED BY DIMENSIONLESS RATIOS; a SETS THE PHYSICAL SCALE.
//
// CBS PHYSICAL IDS (tools/gmsh_to_cbs.py):
//    10 fluid     20 structure solid     21 armour solid ("solid_armour")
//   511 inlet port    520 outlet port    532 heated front face
//   530 every other external surface (adiabatic)
//   Fluid-solid and armour-structure interfaces are conformal (shared nodes).
//   NOTE: the converter currently maps every volume named "solid*" to solid
//   material 1; volume 21 is kept separate so the armour can get its own
//   conductivity once multi-solid materials are added.
//
// USAGE:
//   gmsh ManifoldPanel_CHT.geo -3 -setnumber hScale 2 -o panel.msh
//   gmsh ManifoldPanel_CHT.geo -3 -setnumber nCh 4 -setnumber manifold 1 ...
//==============================================================================

SetFactory("OpenCASCADE");
Geometry.OCCBooleanPreserveNumbering = 1;

//==============================================================================
// 1. PARAMETERS (DIMENSIONLESS RATIOS x CHANNEL SIDE a)
//==============================================================================

DefineConstant[
    a        = {0.010, Name "Geometry/01 channel side a (m)"},
    nCh      = {8,     Name "Geometry/02 number of channels N"},
    pR       = {2.0,   Name "Geometry/03 channel pitch / a"},
    LR       = {30.0,  Name "Geometry/04 channel length / a"},
    wHR      = {2.0,   Name "Geometry/05 header width / a"},
    hHR      = {1.0,   Name "Geometry/06 header height / a"},
    mR       = {1.0,   Name "Geometry/07 solid margin / a"},
    tFrontR  = {0.3,   Name "Geometry/08 front wall / a"},
    tBackR   = {0.7,   Name "Geometry/09 back wall / a"},
    tArmR    = {0.2,   Name "Geometry/10 armour layer / a (0 = none)"},
    manifold = {0,     Name "Geometry/11 manifold: 0 U-type, 1 Z-type"}
];

p      = pR * a;
Lc     = LR * a;
wH     = wHR * a;
hH     = hHR * a;
mS     = mR * a;
tFront = tFrontR * a;
tBack  = tBackR * a;
tArm   = tArmR * a;

// CHANNEL AND HEADER HEIGHTS: THE HEADER IS CENTRED ON THE CHANNEL LAYER AND
// MAY BE TALLER (hH >= a) TO REDUCE ITS PRESSURE LOSS.
zcl  = tBack + hH / 2.0;            // MID-HEIGHT OF THE COOLANT LAYER
zChLo = zcl - a / 2.0;
zChHi = zcl + a / 2.0;
zHLo  = zcl - hH / 2.0;
zHHi  = zcl + hH / 2.0;
Tst   = tBack + hH + tFront;        // STRUCTURE THICKNESS (TO ITS FRONT FACE)
Ttop  = Tst + tArm;                 // HEATED FACE

Wp = nCh * p + 2.0 * mS;            // PANEL WIDTH (y)
Lx = 2.0 * mS + 2.0 * wH + Lc;      // PANEL LENGTH (x)

xHin0  = mS;                        // INLET HEADER
xHin1  = mS + wH;
xHout0 = mS + wH + Lc;              // OUTLET HEADER
xHout1 = xHout0 + wH;

yFirst = -Wp / 2.0 + mS + p / 2.0;  // FIRST CHANNEL CENTRE
yLast  = yFirst + (nCh - 1) * p;    // LAST CHANNEL CENTRE

eps = 1.0e-4 * a;                   // ENTITY-SELECTION TOLERANCE ONLY

If (nCh < 1 || nCh != Round(nCh))
    Error("nCh must be a positive integer (got %g)", nCh);
    Abort;
EndIf
If (pR <= 1.0)
    Error("pitch/a = %g must exceed 1 (channels would touch)", pR);
    Abort;
EndIf
If (hHR < 1.0)
    Error("header height/a = %g must be at least 1", hHR);
    Abort;
EndIf

//==============================================================================
// 2. FLUID: N CHANNELS + TWO HEADERS
//==============================================================================

// CHANNELS RUN INTO THE MIDDLE OF EACH HEADER SO THE UNION IS CLEAN.
vCh[] = {};
For i In {0 : nCh - 1}
    yi = yFirst + i * p;
    vTmp = newv;
    Box(vTmp) = {xHin0 + wH / 2.0, yi - a / 2.0, zChLo, Lc + wH, a, a};
    vCh[] += vTmp;
EndFor

// INLET HEADER: FROM THE PORT FACE y = -Wp/2 TO THE FAR SIDE OF THE LAST CHANNEL.
vHin = newv;
Box(vHin) = {xHin0, -Wp / 2.0, zHLo, wH, (yLast + a / 2.0) + Wp / 2.0, hH};

// OUTLET HEADER: PORT ON THE SAME SIDE (U-TYPE) OR THE OPPOSITE SIDE (Z-TYPE).
vHout = newv;
If (manifold == 0)
    Box(vHout) = {xHout0, -Wp / 2.0, zHLo, wH, (yLast + a / 2.0) + Wp / 2.0, hH};
Else
    Box(vHout) = {xHout0, yFirst - a / 2.0, zHLo, wH, Wp / 2.0 - (yFirst - a / 2.0), hH};
EndIf

vFluidRaw[] = BooleanUnion { Volume{vHin}; Delete; } { Volume{vCh[], vHout}; Delete; };

//==============================================================================
// 3. SOLIDS: STRUCTURE + OPTIONAL ARMOUR, THEN THE CONFORMAL SPLIT
//==============================================================================

vSt = newv;
Box(vSt) = {0.0, -Wp / 2.0, 0.0, Lx, Wp, Tst};
vSolidRaw[] = {vSt};
If (tArm > 0)
    vArm = newv;
    Box(vArm) = {0.0, -Wp / 2.0, Tst, Lx, Wp, tArm};
    vSolidRaw[] += vArm;
EndIf

allVol[] = BooleanFragments { Volume{vSolidRaw[]}; Delete; } { Volume{vFluidRaw[]}; Delete; };

//==============================================================================
// 4. VOLUME IDENTIFICATION
//==============================================================================

// FLUID: EVERYTHING INSIDE THE COOLANT LAYER (SOLIDS SPAN MORE OF z).
// OCC BOUNDING BOXES CARRY A SMALL TOLERANCE, HENCE THE PAD IN z (THE SOLIDS
// EXTEND AT LEAST tBack BELOW AND tFront ABOVE, SO THEY CANNOT FIT).
padz = 0.5 * (tFront < tBack ? tFront : tBack);
vFluid[] = Volume In BoundingBox {
    -eps, -Wp / 2.0 - eps, zHLo - padz,
    Lx + eps, Wp / 2.0 + eps, zHHi + padz };

vArmour[] = {};
If (tArm > 0)
    vArmour[] = Volume In BoundingBox {
        -eps, -Wp / 2.0 - eps, Tst - 0.5 * tArm,
        Lx + eps, Wp / 2.0 + eps, Ttop + eps };
EndIf

vStruct[] = allVol[];
vStruct[] -= vFluid[];
vStruct[] -= vArmour[];

If (#vFluid[] != 1 || #vStruct[] < 1)
    Error("Volume identification failed: fluid %g, structure %g", #vFluid[], #vStruct[]);
    Abort;
EndIf

Physical Volume("fluid", 10) = {vFluid[]};
Physical Volume("solid", 20) = {vStruct[]};
If (tArm > 0)
    Physical Volume("solid_armour", 21) = {vArmour[]};
EndIf

//==============================================================================
// 5. BOUNDARY SURFACES
//==============================================================================

sInlet[] = Surface In BoundingBox {
    xHin0 - eps, -Wp / 2.0 - eps, zHLo - eps,
    xHin1 + eps, -Wp / 2.0 + eps, zHHi + eps };

If (manifold == 0)
    sOutlet[] = Surface In BoundingBox {
        xHout0 - eps, -Wp / 2.0 - eps, zHLo - eps,
        xHout1 + eps, -Wp / 2.0 + eps, zHHi + eps };
Else
    sOutlet[] = Surface In BoundingBox {
        xHout0 - eps, Wp / 2.0 - eps, zHLo - eps,
        xHout1 + eps, Wp / 2.0 + eps, zHHi + eps };
EndIf

sHeated[] = Surface In BoundingBox {
    -eps, -Wp / 2.0 - eps, Ttop - eps,
    Lx + eps, Wp / 2.0 + eps, Ttop + eps };

sExternal[] = Abs(CombinedBoundary { Volume{allVol[]}; });
sAdiabatic[] = sExternal[];
sAdiabatic[] -= sInlet[];
sAdiabatic[] -= sOutlet[];
sAdiabatic[] -= sHeated[];

If (#sInlet[] != 1 || #sOutlet[] != 1 || #sHeated[] != 1)
    Error("Boundary identification failed: inlet %g, outlet %g, heated %g",
          #sInlet[], #sOutlet[], #sHeated[]);
    Abort;
EndIf

Physical Surface("inlet", 511)          = {sInlet[]};
Physical Surface("outlet", 520)         = {sOutlet[]};
Physical Surface("heat_flux", 532)      = {sHeated[]};
Physical Surface("adiabatic_wall", 530) = {sAdiabatic[]};

// FLUID-SOLID INTERFACE (MESH SIZING ONLY).
sInterface[] = Abs(Boundary { Volume{vFluid[]}; });
sInterface[] -= sInlet[];
sInterface[] -= sOutlet[];

//==============================================================================
// 6. MESH SIZES (FLUID AND SOLID GRADED SEPARATELY)
//==============================================================================

DefineConstant[
    hScale = {1.0,   Name "Mesh/0 global size factor (2 = coarser)"},
    hWallR = {0.015, Name "Mesh/1 interface size / a"},
    hCoreR = {0.045, Name "Mesh/2 fluid core size / a"},
    hFrontR= {0.050, Name "Mesh/3 front wall + armour size / a"},
    hSolidR= {0.150, Name "Mesh/4 solid bulk size / a"},
    dWallR = {0.010, Name "Mesh/5 keep interface size up to distance / a"},
    dCoreR = {0.200, Name "Mesh/6 fluid reaches core size at distance / a"},
    dSolidR= {0.800, Name "Mesh/7 solid reaches bulk size at distance / a"}
];

hWall  = hScale * hWallR  * a;
hCore  = hScale * hCoreR  * a;
hFront = hScale * hFrontR * a;
hSolid = hScale * hSolidR * a;

// GRADING DISTANCES SCALE WITH hScale TOO, SO EVERY RESOLUTION GROWS AT THE
// SAME RATE (A COARSE PREVIEW HAS THE SAME ELEMENT QUALITY AS PRODUCTION).
dWall  = hScale * dWallR  * a;
dCore  = hScale * dCoreR  * a;
dSolid = hScale * dSolidR * a;

Field[1] = Distance;
Field[1].SurfacesList = {sInterface[]};
Field[1].Sampling = 60;

Field[2] = Threshold;                 // FLUID: WALL -> CORE
Field[2].InField = 1;
Field[2].SizeMin = hWall;
Field[2].SizeMax = hCore;
Field[2].DistMin = dWall;
Field[2].DistMax = dCore;

Field[3] = Restrict;
Field[3].InField = 2;
Field[3].VolumesList = {vFluid[]};
Field[3].IncludeBoundary = 1;

Field[4] = Threshold;                 // SOLID: WALL -> BULK
Field[4].InField = 1;
Field[4].SizeMin = hWall;
Field[4].SizeMax = hSolid;
Field[4].DistMin = dWall;
Field[4].DistMax = dSolid;

Field[5] = Restrict;
Field[5].InField = 4;
Field[5].VolumesList = {vStruct[], vArmour[]};
Field[5].IncludeBoundary = 1;

Field[6] = Box;                       // FRONT WALL + ARMOUR: STEEPEST GRADIENTS
Field[6].VIn = hFront;
Field[6].VOut = hSolid;
Field[6].XMin = -eps;
Field[6].XMax = Lx + eps;
Field[6].YMin = -Wp / 2.0 - eps;
Field[6].YMax =  Wp / 2.0 + eps;
Field[6].ZMin = zHHi;
Field[6].ZMax = Ttop + eps;
Field[6].Thickness = tFront;

Field[7] = Min;
Field[7].FieldsList = {3, 5, 6};
Background Field = 7;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

//==============================================================================
// 7. TRUE 3-D LINEAR TETRAHEDRAL MESH
//==============================================================================

DefineConstant[ algo3D = {10, Name "Mesh/8 3-D algorithm (10 HXT, 1 Delaunay)"} ];

Mesh.Algorithm = 6;
Mesh.Algorithm3D = algo3D;
Mesh.ElementOrder = 1;
Mesh.Smoothing = 10;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.OptimizeThreshold = 0.35;
Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;
Mesh.SaveAll = 0;

//==============================================================================
// 8. AUDIT
//==============================================================================

Printf("==============================================================");
Printf("MANIFOLD-FED COOLED PANEL, manifold type %g (0 = U, 1 = Z)", manifold);
Printf("channels N, a, pitch    = %g, %.4e m, %.4e m", nCh, a, p);
Printf("channel length          = %.4e m (%g a)", Lc, LR);
Printf("headers w x h           = %.4e x %.4e m", wH, hH);
Printf("panel Lx x Wp x Ttop    = %.4e x %.4e x %.4e m", Lx, Wp, Ttop);
Printf("back / front / armour   = %.4e / %.4e / %.4e m", tBack, tFront, tArm);
Printf("volumes fluid/struct/arm = %g / %g / %g", #vFluid[], #vStruct[], #vArmour[]);
Printf("faces: inlet %g, outlet %g, heated %g, adiabatic %g, interface %g",
       #sInlet[], #sOutlet[], #sHeated[], #sAdiabatic[], #sInterface[]);
Printf("sizes: wall %.3e, core %.3e, front %.3e, solid %.3e m",
       hWall, hCore, hFront, hSolid);
Printf("==============================================================");
