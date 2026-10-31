//==============================================================================
// BlanketSlice_CHT.geo
//
// BLANKET SLICE FOR CBS3D++_SI CONJUGATE HEAT TRANSFER:
// A TAPERED SOLID PLATE WITH A FIRST WALL AND 4 NESTED U-SHAPED COOLANT
// CHANNELS, PORTS ON THE BACK FACE, HEATED ON THE FRONT (FIRST-WALL) FACE.
//
// DIMENSIONS FOLLOW THE DESIGN DRAWING (mm, x FROM THE LEFT BACK CORNER):
//
//   PLATE      back width 350, height 280, sides at 85 deg, front fillet R6
//   CHANNEL 1  top edge 270 (first wall 10), width 15, leg at x = 20, R10
//   CHANNEL 2  top edge 220,                 width 10, leg at x = 65, R5
//   CHANNEL 3  top edge 150,                 width 10, leg at x = 105, R5
//   CHANNEL 4  top edge 100,                 width  5, leg at x = 135, R3.5
//
//   "top edge" = OUTER (UPPER) EDGE OF THE CHANNEL'S FRONT RUN, FROM y = 0.
//   "leg at x" = OUTER EDGE OF THE LEFT LEG ON THE BACK FACE y = 0.
//   Every channel leg is parallel to the plate side (85 deg). Both edges of a
//   channel corner carry the SAME fillet radius, as drawn. Channel widths are
//   true (perpendicular) widths. The plate is symmetric about its centre line.
//
//                 R6  ______ FIRST WALL (heated) ______  R6
//                    /  R10 ____________________ R10    \
//                   /  /  ________________________ \  \   \
//                  /  /  /  R5 ______________ R5  \  \  \  \
//                 /  /  /  /  ______________   \  \  \  \  \
//                /  /  /  /  /  R3.5 ___ R3.5 \  \  \  \  \  \
//               /  /  /  /  /  /         \  \  \  \  \  \  \  \
//              /__|__|__|__|__|___SOLID___|__|__|__|__|__|__|__\  y = 0
//                 inlets (x < centre)      outlets (x > centre)
//
// THROUGH-THICKNESS (z): EACH CHANNEL IS SQUARE (DEPTH = WIDTH) AND CENTRED
// AT z = Tz/2; Tz = w1 + 2 tz.
//
// CBS PHYSICAL IDS (tools/gmsh_to_cbs.py):
//    10 fluid   20 solid   21 armour (only if tArm > 0)
//   511 inlets (4 left leg ends)   520 outlets (4 right leg ends)
//   532 heated first-wall face (front face + both R6 fillets)
//   530 every other external face (adiabatic)
//
// COORDINATES IN THIS FILE: x = 0 ON THE CENTRE LINE (DRAWING x MINUS 175 mm),
// y = 0 ON THE BACK FACE, LENGTHS IN METRES.
//
// USAGE:
//   gmsh BlanketSlice_CHT.geo                                       geometry
//   gmsh BlanketSlice_CHT.geo -3 -setnumber hScale 4 -o slice.msh   coarse mesh
//==============================================================================

SetFactory("OpenCASCADE");
Geometry.OCCBooleanPreserveNumbering = 1;

//==============================================================================
// 1. PARAMETERS (METRES)
//==============================================================================

DefineConstant[
    Wb    = {0.350,  Name "Plate/1 back width (m)"},
    H     = {0.280,  Name "Plate/2 height, back face to front (m)"},
    angD  = {85.0,   Name "Plate/3 side angle to the back face (deg)"},
    rP    = {0.006,  Name "Plate/4 front fillet radius (m)"},

    yT1   = {0.270,  Name "Channel 1/1 top edge height (m)"},
    w1    = {0.015,  Name "Channel 1/2 width (m)"},
    x1    = {0.020,  Name "Channel 1/3 leg outer edge, x from left corner (m)"},
    r1    = {0.010,  Name "Channel 1/4 fillet radius (m)"},

    yT2   = {0.220,  Name "Channel 2/1 top edge height (m)"},
    w2    = {0.010,  Name "Channel 2/2 width (m)"},
    x2    = {0.065,  Name "Channel 2/3 leg outer edge, x from left corner (m)"},
    r2    = {0.005,  Name "Channel 2/4 fillet radius (m)"},

    yT3   = {0.150,  Name "Channel 3/1 top edge height (m)"},
    w3    = {0.010,  Name "Channel 3/2 width (m)"},
    x3    = {0.105,  Name "Channel 3/3 leg outer edge, x from left corner (m)"},
    r3    = {0.005,  Name "Channel 3/4 fillet radius (m)"},

    yT4   = {0.100,  Name "Channel 4/1 top edge height (m)"},
    w4    = {0.005,  Name "Channel 4/2 width (m)"},
    x4    = {0.135,  Name "Channel 4/3 leg outer edge, x from left corner (m)"},
    r4    = {0.0035, Name "Channel 4/4 fillet radius (m)"},

    tzR   = {0.5,    Name "Thickness/1 z wall each side of channel 1 / w1"},
    squareCh = {1,   Name "Thickness/2 square channels (1) or all depth w1 (0)"},
    tArm  = {0.0,    Name "Thickness/3 armour on the first wall (m), 0 = none"}
];

nCh = 4;
yT[] = {yT1, yT2, yT3, yT4};
wC[] = {w1,  w2,  w3,  w4};
xC[] = {x1,  x2,  x3,  x4};
rF[] = {r1,  r2,  r3,  r4};

ang = angD * Pi / 180.0;
tz  = tzR * w1;
Tz  = w1 + 2.0 * tz;                     // SLICE THICKNESS
zMid = Tz / 2.0;
xL0 = -Wb / 2.0;                         // LEFT BACK CORNER (CENTRED COORDS)

// LEFT SIDE DIRECTION (UP AND INWARD) AND INWARD NORMAL.
tx = Cos(ang);  ty = Sin(ang);
nx = Sin(ang);  ny = -Cos(ang);

eps = 1.0e-6;                            // ENTITY-SELECTION TOLERANCE ONLY

//==============================================================================
// 2. MACRO: REGION BOUNDED BY y = 0, TWO SIDES PARALLEL TO THE PLATE SIDES
//    (LEFT SIDE THROUGH (xF, 0), MIRRORED ABOUT x = 0) AND A FRONT LINE
//    y = yTop, WITH FILLETS OF RADIUS rr AT THE TWO FRONT CORNERS.
//    FLAT FACE AT z = zB, EXTRUDED BY hZ.
//    INPUT : xF, yTop, rr, zB, hZ        OUTPUT: vReg (VOLUME TAG)
//==============================================================================

Macro FilletedRegion
    // FILLET CENTRE: rr INSIDE THE SIDE LINE AND rr BELOW THE FRONT LINE.
    cyR = yTop - rr;
    sR  = (cyR - rr * ny) / ty;
    cxR = xF + rr * nx + sR * tx;

    tsxR = cxR - rr * nx;                // TANGENT POINT ON THE SIDE
    tsyR = cyR - rr * ny;
    tfxR = cxR;                          // TANGENT POINT ON THE FRONT
    tfyR = yTop;

    If (cxR >= 0.0 || tsyR <= 0.0 || rr <= 0.0)
        Error("Region does not fit: foot x %g, top %g, radius %g", xF, yTop, rr);
        Abort;
    EndIf

    pA = newp; Point(pA) = { xF,   0.0,  zB};
    pB = newp; Point(pB) = { tsxR, tsyR, zB};
    pC = newp; Point(pC) = { cxR,  cyR,  zB};
    pD = newp; Point(pD) = { tfxR, tfyR, zB};
    pE = newp; Point(pE) = {-tfxR, tfyR, zB};
    pF = newp; Point(pF) = {-cxR,  cyR,  zB};
    pG = newp; Point(pG) = {-tsxR, tsyR, zB};
    pH = newp; Point(pH) = {-xF,   0.0,  zB};

    lAB = newl; Line(lAB)   = {pA, pB};
    lBD = newl; Circle(lBD) = {pB, pC, pD};
    lDE = newl; Line(lDE)   = {pD, pE};
    lEG = newl; Circle(lEG) = {pE, pF, pG};
    lGH = newl; Line(lGH)   = {pG, pH};
    lHA = newl; Line(lHA)   = {pH, pA};

    clR = newcl; Curve Loop(clR) = {lAB, lBD, lDE, lEG, lGH, lHA};
    sRg = news;  Plane Surface(sRg) = {clR};
    exR[] = Extrude {0, 0, hZ} { Surface{sRg}; };
    vReg = exR[1];
Return

//==============================================================================
// 3. PLATE
//==============================================================================

xF = xL0; yTop = H; rr = rP; zB = 0.0; hZ = Tz;
Call FilletedRegion;
vPlate = vReg;

//==============================================================================
// 4. CHANNELS: OUTER REGION MINUS INNER REGION (SAME FILLET RADIUS ON BOTH
//    EDGES, AS DRAWN). THE INNER REGION IS THE OUTER ONE MOVED IN BY THE
//    CHANNEL WIDTH: SIDES BY w (PERPENDICULAR), FRONT DOWN BY w.
//==============================================================================

vChan[] = {};
For i In {0 : nCh - 1}
    w  = wC[i];
    dZ = (squareCh == 1) ? w : w1;       // CHANNEL DEPTH IN z

    xF = xL0 + xC[i]; yTop = yT[i]; rr = rF[i]; zB = zMid - dZ / 2.0; hZ = dZ;
    Call FilletedRegion;
    vOuter = vReg;

    xF = xL0 + xC[i] + w / nx; yTop = yT[i] - w; rr = rF[i]; zB = zMid - dZ / 2.0; hZ = dZ;
    Call FilletedRegion;
    vInner = vReg;

    vBand[] = BooleanDifference { Volume{vOuter}; Delete; } { Volume{vInner}; Delete; };
    vChan[] += vBand[];

    Printf("channel %g: width %.4e m, depth %.4e m, top edges %.4e / %.4e m, leg x %.4e m, fillet %.4e m",
           i + 1, w, dZ, yT[i], yT[i] - w, xC[i], rF[i]);
EndFor

//==============================================================================
// 5. OPTIONAL ARMOUR STRIP AND CONFORMAL SPLIT
//==============================================================================

vSolidParts[] = {vPlate};
If (tArm > 0)
    vPlateCopy[] = Translate {0, 0, 0} { Duplicata { Volume{vPlate}; } };
    bA = newv; Box(bA) = {-Wb, H - tArm, -Tz, 2.0 * Wb, 2.0 * tArm, 3.0 * Tz};
    bB = newv; Box(bB) = {-Wb, H - tArm, -Tz, 2.0 * Wb, 2.0 * tArm, 3.0 * Tz};
    vArmRaw[] = BooleanIntersection { Volume{vPlateCopy[]}; Delete; } { Volume{bA}; Delete; };
    vStRaw[]  = BooleanDifference   { Volume{vPlate};       Delete; } { Volume{bB}; Delete; };
    vSolidParts[] = {vStRaw[], vArmRaw[]};
EndIf

allVol[] = BooleanFragments { Volume{vSolidParts[]}; Delete; } { Volume{vChan[]}; Delete; };

//==============================================================================
// 6. VOLUME IDENTIFICATION
//==============================================================================

padz = 0.5 * tz;                         // CHANNELS STAY INSIDE THE w1 LAYER
vFluid[] = Volume In BoundingBox {
    -Wb, -eps, zMid - w1 / 2.0 - padz,
     Wb,  H + eps, zMid + w1 / 2.0 + padz };

vArmour[] = {};
If (tArm > 0)
    vArmour[] = Volume In BoundingBox { -Wb, H - 1.5 * tArm, -eps, Wb, H + eps, Tz + eps };
EndIf

vSolid[] = allVol[];
vSolid[] -= vFluid[];
vSolid[] -= vArmour[];

If (#vFluid[] != nCh || #vSolid[] < 1)
    Error("Volume identification failed: %g fluid (expected %g), %g solid",
          #vFluid[], nCh, #vSolid[]);
    Abort;
EndIf

Physical Volume("fluid", 10) = {vFluid[]};
Physical Volume("solid", 20) = {vSolid[]};
If (tArm > 0)
    Physical Volume("solid_armour", 21) = {vArmour[]};
EndIf

//==============================================================================
// 7. BOUNDARY SURFACES
//==============================================================================

sInlet[] = Surface In BoundingBox {
    -Wb, -eps, zMid - w1 / 2.0 - eps,
    0.0,  eps, zMid + w1 / 2.0 + eps };
sOutlet[] = Surface In BoundingBox {
    0.0, -eps, zMid - w1 / 2.0 - eps,
     Wb,  eps, zMid + w1 / 2.0 + eps };

// HEATED FIRST WALL: EXTERNAL FACES THAT REACH y = H (FRONT FACE AND THE
// R6 FILLETS), EXCLUDING THE TWO z FACES.
sExternal[] = Abs(CombinedBoundary { Volume{allVol[]}; });
sHeated[] = {};
For k In {0 : #sExternal[] - 1}
    bbS[] = BoundingBox Surface{ sExternal[k] };
    isZFace = (bbS[5] - bbS[2] < eps);
    If (bbS[4] >= H - eps && isZFace == 0)
        sHeated[] += sExternal[k];
    EndIf
EndFor

sAdiabatic[] = sExternal[];
sAdiabatic[] -= sInlet[];
sAdiabatic[] -= sOutlet[];
sAdiabatic[] -= sHeated[];

If (#sInlet[] != nCh || #sOutlet[] != nCh || #sHeated[] < 3)
    Error("Boundary identification failed: inlets %g, outlets %g (expected %g), heated %g",
          #sInlet[], #sOutlet[], nCh, #sHeated[]);
    Abort;
EndIf

Physical Surface("inlet", 511)          = {sInlet[]};
Physical Surface("outlet", 520)         = {sOutlet[]};
Physical Surface("heat_flux", 532)      = {sHeated[]};
Physical Surface("adiabatic_wall", 530) = {sAdiabatic[]};

sInterface[] = Abs(Boundary { Volume{vFluid[]}; });
sInterface[] -= sInlet[];
sInterface[] -= sOutlet[];

//==============================================================================
// 8. MESH SIZES (RATIOS OF THE SMALLEST CHANNEL WIDTH w4). FLUID AND SOLID
//    ARE GRADED SEPARATELY; SIZES AND GRADING DISTANCES SCALE WITH hScale.
//==============================================================================

DefineConstant[
    hScale  = {1.0,  Name "Mesh/0 global size factor (2 = coarser)"},
    hWallR  = {0.03, Name "Mesh/1 interface size / w4"},
    hCoreR  = {0.08, Name "Mesh/2 fluid core size / w4"},
    hFrontR = {0.12, Name "Mesh/3 first wall size / w4"},
    hSolidR = {0.50, Name "Mesh/4 solid bulk size / w4"},
    dWallR  = {0.02, Name "Mesh/5 keep interface size up to distance / w4"},
    dCoreR  = {0.40, Name "Mesh/6 fluid reaches core size at distance / w4"},
    dSolidR = {2.00, Name "Mesh/7 solid reaches bulk size at distance / w4"}
];

hWall  = hScale * hWallR  * w4;
hCore  = hScale * hCoreR  * w4;
hFront = hScale * hFrontR * w4;
hSolid = hScale * hSolidR * w4;
dWall  = hScale * dWallR  * w4;
dCoreS = hScale * dCoreR  * w4;
dSolid = hScale * dSolidR * w4;
tFW    = H - yT1;                        // FIRST-WALL THICKNESS

Field[1] = Distance;
Field[1].SurfacesList = {sInterface[]};
Field[1].Sampling = 60;

Field[2] = Threshold;                    // FLUID: WALL -> CORE
Field[2].InField = 1;
Field[2].SizeMin = hWall;
Field[2].SizeMax = hCore;
Field[2].DistMin = dWall;
Field[2].DistMax = dCoreS;

Field[3] = Restrict;
Field[3].InField = 2;
Field[3].VolumesList = {vFluid[]};
Field[3].IncludeBoundary = 1;

Field[4] = Threshold;                    // SOLID: WALL -> BULK
Field[4].InField = 1;
Field[4].SizeMin = hWall;
Field[4].SizeMax = hSolid;
Field[4].DistMin = dWall;
Field[4].DistMax = dSolid;

Field[5] = Restrict;
Field[5].InField = 4;
Field[5].VolumesList = {vSolid[], vArmour[]};
Field[5].IncludeBoundary = 1;

Field[6] = Distance;                     // FIRST WALL: FINE THROUGH ITS THICKNESS
Field[6].SurfacesList = {sHeated[]};
Field[6].Sampling = 60;

Field[7] = Threshold;
Field[7].InField = 6;
Field[7].SizeMin = hFront;
Field[7].SizeMax = hSolid;
Field[7].DistMin = tFW;
Field[7].DistMax = tFW + dSolid;

Field[8] = Restrict;
Field[8].InField = 7;
Field[8].VolumesList = {vSolid[], vArmour[]};
Field[8].IncludeBoundary = 1;

Field[9] = Min;
Field[9].FieldsList = {3, 5, 8};
Background Field = 9;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

//==============================================================================
// 9. TRUE 3-D LINEAR TETRAHEDRAL MESH
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
// 10. AUDIT
//==============================================================================

Printf("==============================================================");
Printf("BLANKET SLICE: PLATE %.4e (back) x %.4e (height) m, sides %g deg, R %.4e m",
       Wb, H, angD, rP);
Printf("front width = %.4e m, first wall = %.4e m, slice thickness Tz = %.4e m",
       Wb - 2.0 * H / Tan(ang), tFW, Tz);
Printf("volumes fluid / solid / armour = %g / %g / %g", #vFluid[], #vSolid[], #vArmour[]);
Printf("faces: inlets %g, outlets %g, heated %g, adiabatic %g, interface %g",
       #sInlet[], #sOutlet[], #sHeated[], #sAdiabatic[], #sInterface[]);
Printf("sizes: wall %.3e, core %.3e, first wall %.3e, solid %.3e m",
       hWall, hCore, hFront, hSolid);
Printf("==============================================================");
