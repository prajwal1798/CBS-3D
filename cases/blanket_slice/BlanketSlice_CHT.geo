//==============================================================================
// BlanketSlice_CHT.geo
//
// BLANKET SLICE: A TAPERED COOLING PLATE WITH 4 NESTED U-SHAPED SQUARE
// CHANNELS UNDER A FIRST WALL WITH ARMOUR, FOR CBS3D++_SI CONJUGATE HEAT
// TRANSFER.
//
// IN-PLANE GEOMETRY (x = WIDTH, y = DEPTH FROM THE BACK FACE; THE PLATE IS A
// SLICE OF THICKNESS Tz IN z):
//
//                 armour + first wall   <== HEAT FLUX (front + fillets)
//                   ______________________
//                .'  __________________  '.          channel 1: largest, next
//              .'  .'  ______________  '.  '.        to the first wall
//            .'  .'  .'  __________  '.  '.  '.      channels 2-4: smaller,
//           /  .'  .'  .'          '.  '.  '.  \     equal gaps between them
//          /  /  /  /    solid core   \  \  \  \
//         /  /  /  /    (cut out)      \  \  \  \
//        /__|__|__|____________________|__|__|__\    y = 0 (back face)
//           ^  ^  ^  ^              v  v  v  v
//           inlets (x < 0)          outlets (x > 0)
//
// UNIFORM-OFFSET RULE: THE PLATE OUTLINE, BOTH EDGES OF EVERY CHANNEL AND THE
// CORE CUT-OUT ARE ALL OFFSETS OF ONE ROUNDED TAPERED OUTLINE. AT EACH FRONT
// CORNER EVERY ARC SHARES ONE CENTRE, SO
//   - CORNERS ARE PARALLEL (CONCENTRIC) THROUGHOUT,
//   - EACH CHANNEL KEEPS A CONSTANT WIDTH THROUGH ITS BEND,
//   - THE GAP BETWEEN NEIGHBOURING CHANNELS IS CONSTANT EVERYWHERE.
// AN ARC AT OFFSET d FROM THE OUTLINE HAS RADIUS R0 - d, SO R0 MUST EXCEED
// THE DEEPEST OFFSET (CHECKED BELOW).
//
// THROUGH-THICKNESS (z): EACH CHANNEL IS SQUARE (w x w), CENTRED AT z = Tz/2.
// Tz = w1 + 2 tz.
//
// CBS PHYSICAL IDS (tools/gmsh_to_cbs.py):
//    10 fluid     20 structure     21 armour ("solid_armour")
//   511 inlets (left leg ends)    520 outlets (right leg ends)
//   532 heated front (front face + both front fillets)
//   530 all other external faces (adiabatic)
//
// USAGE (EVERY PARAMETER CAN BE SET FROM THE COMMAND LINE):
//   gmsh BlanketSlice_CHT.geo                                      geometry only
//   gmsh BlanketSlice_CHT.geo -3 -setnumber hScale 5 -o slice.msh  coarse mesh
//==============================================================================

SetFactory("OpenCASCADE");
Geometry.OCCBooleanPreserveNumbering = 1;

//==============================================================================
// 1. PARAMETERS (LENGTHS IN METRES; RATIOS ARE x a = CHANNEL 1 WIDTH)
//==============================================================================

DefineConstant[
    a     = {0.012, Name "Geometry/01 channel 1 width a (m)"},
    Wb    = {0.350, Name "Geometry/02 back width at y = 0 (m)"},
    Wf    = {0.270, Name "Geometry/03 front width at y = H (m)"},
    H     = {0.400, Name "Geometry/04 depth H, back face to front (m)"},
    R0    = {0.100, Name "Geometry/05 front fillet radius of the plate (m)"},
    w2R   = {0.80,  Name "Geometry/06 channel 2 width / a"},
    w3R   = {0.65,  Name "Geometry/07 channel 3 width / a"},
    w4R   = {0.50,  Name "Geometry/08 channel 4 width / a"},
    gR    = {1.00,  Name "Geometry/09 gap between channels / a"},
    tArmR = {0.20,  Name "Geometry/10 armour thickness / a (0 = none)"},
    tFWR  = {0.30,  Name "Geometry/11 first wall (armour to channel 1) / a"},
    tzR   = {0.50,  Name "Geometry/12 z wall each side of channel 1 / a"},
    coreCut = {1,   Name "Geometry/13 remove solid inside channel 4 (1 = yes)"},
    gCoreR  = {1.00,Name "Geometry/14 wall left inside channel 4 / a"}
];

nCh  = 4;
wCh[] = {a, w2R * a, w3R * a, w4R * a};    // CHANNEL WIDTHS (SQUARE SECTIONS)
g    = gR * a;
tArm = tArmR * a;
tFW  = tFWR * a;
tz   = tzR * a;
Tz   = a + 2.0 * tz;                       // SLICE THICKNESS
zMid = Tz / 2.0;
gCore = gCoreR * a;

// CENTRE-LINE OFFSET OF EACH CHANNEL FROM THE PLATE OUTLINE.
dCh[] = {};
dEdge = tArm + tFW;                        // OUTER EDGE OF CHANNEL 1
For i In {0 : nCh - 1}
    dCh[] += dEdge + wCh[i] / 2.0;
    dEdge = dEdge + wCh[i] + g;
EndFor
dInner = dCh[nCh - 1] + wCh[nCh - 1] / 2.0;    // INNER EDGE OF CHANNEL 4
dCore  = dInner + gCore;                       // CORE CUT-OUT BOUNDARY

// LEFT SIDE: FROM (-Wb/2, 0) TO (-Wf/2, H). UNIT TANGENT t, INWARD NORMAL n.
dxs = (Wb - Wf) / 2.0;
Ls  = Sqrt(dxs * dxs + H * H);
tx  = dxs / Ls;
ty  = H / Ls;
nx  = H / Ls;
ny  = -dxs / Ls;

eps = 1.0e-4 * a;                          // ENTITY-SELECTION TOLERANCE ONLY

dDeepest = (coreCut == 1) ? dCore : dInner;
If (Wf >= Wb || Wf <= 0)
    Error("Need 0 < Wf < Wb (got Wf = %g, Wb = %g)", Wf, Wb);
    Abort;
EndIf
If (R0 - dDeepest < 0.25 * a)
    Error("R0 = %g is too small: the deepest offset is %g, so its arc radius would be %g (< a/4). Increase R0 or reduce gaps/widths.",
          R0, dDeepest, R0 - dDeepest);
    Abort;
EndIf

//==============================================================================
// 2. MACRO: OFFSET REGION. THE PLATE OUTLINE OFFSET INWARD BY dOff, WITH THE
//    CONCENTRIC CORNER RADIUS R0 - dOff, AS A FLAT SURFACE AT z = zB EXTRUDED
//    BY hZ. THE REGION IS OPEN TO THE BACK FACE y = 0.
//    INPUT : dOff, zB, hZ        OUTPUT: vTrap (VOLUME TAG)
//==============================================================================

Macro OffsetRegion
    rC = R0 - dOff;

    // OFFSET SIDE LINE: POINT (-Wb/2, 0) + dOff n, DIRECTION t; ITS FOOT ON y = 0.
    q0x = -Wb / 2.0 + dOff * nx;
    q0y = dOff * ny;
    xBot = q0x + (-q0y / ty) * tx;

    // COMMON CORNER CENTRE: DISTANCE R0 FROM BOTH THE SIDE AND THE FRONT.
    k0x = -Wb / 2.0 + R0 * nx;
    k0y = R0 * ny;
    cyC = H - R0;
    cxC = k0x + ((cyC - k0y) / ty) * tx;

    // TANGENT POINTS ON THE OFFSET SIDE AND ON THE OFFSET FRONT.
    tsx = cxC - rC * nx;
    tsy = cyC - rC * ny;
    tfx = cxC;
    tfy = H - dOff;

    If (cxC >= 0.0 || tsy <= 0.0 || rC <= 0.0)
        Error("Offset region does not fit: offset %g, radius %g (centre x %g, side tangent y %g)",
              dOff, rC, cxC, tsy);
        Abort;
    EndIf

    pA = newp; Point(pA) = { xBot, 0.0, zB};
    pB = newp; Point(pB) = { tsx,  tsy, zB};
    pC = newp; Point(pC) = { cxC,  cyC, zB};
    pD = newp; Point(pD) = { tfx,  tfy, zB};
    pE = newp; Point(pE) = {-tfx,  tfy, zB};
    pF = newp; Point(pF) = {-cxC,  cyC, zB};
    pG = newp; Point(pG) = {-tsx,  tsy, zB};
    pH = newp; Point(pH) = {-xBot, 0.0, zB};

    lAB = newl; Line(lAB)   = {pA, pB};
    lBD = newl; Circle(lBD) = {pB, pC, pD};
    lDE = newl; Line(lDE)   = {pD, pE};
    lEG = newl; Circle(lEG) = {pE, pF, pG};
    lGH = newl; Line(lGH)   = {pG, pH};
    lHA = newl; Line(lHA)   = {pH, pA};

    clT = newcl; Curve Loop(clT) = {lAB, lBD, lDE, lEG, lGH, lHA};
    sT  = news;  Plane Surface(sT) = {clT};
    exT[] = Extrude {0, 0, hZ} { Surface{sT}; };
    vTrap = exT[1];
Return

//==============================================================================
// 3. PLATE, WITH THE OPTIONAL CORE CUT-OUT
//==============================================================================

dOff = 0.0; zB = 0.0; hZ = Tz;
Call OffsetRegion;
vPlate = vTrap;

If (coreCut == 1)
    dOff = dCore; zB = -a; hZ = Tz + 2.0 * a;
    Call OffsetRegion;
    vPlateCut[] = BooleanDifference { Volume{vPlate}; Delete; } { Volume{vTrap}; Delete; };
    vPlate = vPlateCut[0];
EndIf

//==============================================================================
// 4. CHANNELS: BAND BETWEEN TWO CONCENTRIC OFFSETS, SQUARE SECTION w x w
//==============================================================================

vChan[] = {};
For i In {0 : nCh - 1}
    w  = wCh[i];
    dC = dCh[i];

    dOff = dC - w / 2.0; zB = zMid - w / 2.0; hZ = w;
    Call OffsetRegion;
    vOuter = vTrap;

    dOff = dC + w / 2.0; zB = zMid - w / 2.0; hZ = w;
    Call OffsetRegion;
    vInner = vTrap;

    vBand[] = BooleanDifference { Volume{vOuter}; Delete; } { Volume{vInner}; Delete; };
    vChan[] += vBand[];

    Printf("channel %g: %.4e x %.4e m, centre-line offset %.4e m, bend radius %.4e m (inner %.4e, outer %.4e)",
           i + 1, w, w, dC, R0 - dC, R0 - dC - w / 2.0, R0 - dC + w / 2.0);
EndFor

//==============================================================================
// 5. ARMOUR STRIP (y >= H - tArm) AND CONFORMAL SPLIT
//==============================================================================

vSolidParts[] = {vPlate};
If (tArm > 0)
    vPlateCopy[] = Translate {0, 0, 0} { Duplicata { Volume{vPlate}; } };
    bA = newv;
    Box(bA) = {-Wb, H - tArm, -a, 2.0 * Wb, 2.0 * tArm, Tz + 2.0 * a};
    bB = newv;
    Box(bB) = {-Wb, H - tArm, -a, 2.0 * Wb, 2.0 * tArm, Tz + 2.0 * a};
    vArmRaw[] = BooleanIntersection { Volume{vPlateCopy[]}; Delete; } { Volume{bA}; Delete; };
    vStRaw[]  = BooleanDifference   { Volume{vPlate};       Delete; } { Volume{bB}; Delete; };
    vSolidParts[] = {vStRaw[], vArmRaw[]};
EndIf

allVol[] = BooleanFragments { Volume{vSolidParts[]}; Delete; } { Volume{vChan[]}; Delete; };

//==============================================================================
// 6. VOLUME IDENTIFICATION
//==============================================================================

// FLUID: THE CHANNELS LIE INSIDE THE z-RANGE OF CHANNEL 1; THE SOLIDS SPAN Tz.
padz = 0.5 * tz;
vFluid[] = Volume In BoundingBox {
    -Wb, -eps, zMid - a / 2.0 - padz,
     Wb,  H + eps, zMid + a / 2.0 + padz };

vArmour[] = {};
If (tArm > 0)
    vArmour[] = Volume In BoundingBox {
        -Wb, H - 1.5 * tArm, -eps,
         Wb, H + eps, Tz + eps };
EndIf

vStruct[] = allVol[];
vStruct[] -= vFluid[];
vStruct[] -= vArmour[];

If (#vFluid[] != nCh || #vStruct[] < 1)
    Error("Volume identification failed: %g fluid (expected %g), %g structure",
          #vFluid[], nCh, #vStruct[]);
    Abort;
EndIf

Physical Volume("fluid", 10) = {vFluid[]};
Physical Volume("solid", 20) = {vStruct[]};
If (tArm > 0)
    Physical Volume("solid_armour", 21) = {vArmour[]};
EndIf

//==============================================================================
// 7. BOUNDARY SURFACES
//==============================================================================

// PORTS: CHANNEL END FACES ON y = 0 (INSIDE THE CHANNEL LAYER).
sInlet[] = Surface In BoundingBox {
    -Wb, -eps, zMid - a / 2.0 - eps,
    0.0,  eps, zMid + a / 2.0 + eps };
sOutlet[] = Surface In BoundingBox {
    0.0, -eps, zMid - a / 2.0 - eps,
     Wb,  eps, zMid + a / 2.0 + eps };

// HEATED FRONT: FRONT FACE + BOTH FRONT FILLETS (THE CURVED FACES START AT
// THE SIDE TANGENT POINT, BELOW y = H - R0), EXCLUDING THE z FACES.
yHeat = H - R0 - R0 * Fabs(ny) - eps;
sFront[] = Surface In BoundingBox {
    -Wb, yHeat, -eps,
     Wb, H + eps, Tz + eps };
sZ0[] = Surface In BoundingBox { -Wb, -eps, -eps,      Wb, H + eps, eps };
sZ1[] = Surface In BoundingBox { -Wb, -eps, Tz - eps,  Wb, H + eps, Tz + eps };
sCand[] = sFront[];
sCand[] -= sZ0[];
sCand[] -= sZ1[];

// KEEP ONLY EXTERNAL FACES THAT REACH THE FRONT LAYER (y >= H - tArm): THIS
// DROPS CHANNEL WALLS, THE ARMOUR CUT AND THE FACES OF THE CORE CUT-OUT.
sExternal[] = Abs(CombinedBoundary { Volume{allVol[]}; });
sHeated[] = {};
For k In {0 : #sCand[] - 1}
    For m In {0 : #sExternal[] - 1}
        If (sCand[k] == sExternal[m])
            bbS[] = BoundingBox Surface{ sCand[k] };
            If (bbS[4] >= H - tArm - 1.0e-3 * a)
                sHeated[] += sCand[k];
            EndIf
        EndIf
    EndFor
EndFor

sAdiabatic[] = sExternal[];
sAdiabatic[] -= sInlet[];
sAdiabatic[] -= sOutlet[];
sAdiabatic[] -= sHeated[];

If (#sInlet[] != nCh || #sOutlet[] != nCh || #sHeated[] != 5)
    Error("Boundary identification failed: inlets %g, outlets %g (expected %g), heated %g",
          #sInlet[], #sOutlet[], nCh, #sHeated[]);
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
// 8. MESH SIZES: FLUID AND SOLID GRADED SEPARATELY. SIZES ARE RATIOS OF THE
//    SMALLEST CHANNEL WIDTH; SIZES AND GRADING DISTANCES SCALE WITH hScale.
//==============================================================================

DefineConstant[
    hScale  = {1.0,   Name "Mesh/0 global size factor (2 = coarser)"},
    hWallR  = {0.025, Name "Mesh/1 interface size / smallest channel"},
    hCoreR  = {0.070, Name "Mesh/2 fluid core size / smallest channel"},
    hFrontR = {0.080, Name "Mesh/3 first wall + armour size / smallest channel"},
    hSolidR = {0.400, Name "Mesh/4 solid bulk size / smallest channel"},
    dWallR  = {0.015, Name "Mesh/5 keep interface size up to distance / smallest"},
    dCoreR  = {0.300, Name "Mesh/6 fluid reaches core size at distance / smallest"},
    dSolidR = {1.500, Name "Mesh/7 solid reaches bulk size at distance / smallest"}
];

wMin   = wCh[nCh - 1];
hWall  = hScale * hWallR  * wMin;
hCore  = hScale * hCoreR  * wMin;
hFront = hScale * hFrontR * wMin;
hSolid = hScale * hSolidR * wMin;
dWall  = hScale * dWallR  * wMin;
dCoreS = hScale * dCoreR  * wMin;
dSolid = hScale * dSolidR * wMin;

// F1-F3: FLUID, WALL -> CORE.
Field[1] = Distance;
Field[1].SurfacesList = {sInterface[]};
Field[1].Sampling = 60;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = hWall;
Field[2].SizeMax = hCore;
Field[2].DistMin = dWall;
Field[2].DistMax = dCoreS;

Field[3] = Restrict;
Field[3].InField = 2;
Field[3].VolumesList = {vFluid[]};
Field[3].IncludeBoundary = 1;

// F4-F5: SOLIDS, WALL -> BULK.
Field[4] = Threshold;
Field[4].InField = 1;
Field[4].SizeMin = hWall;
Field[4].SizeMax = hSolid;
Field[4].DistMin = dWall;
Field[4].DistMax = dSolid;

Field[5] = Restrict;
Field[5].InField = 4;
Field[5].VolumesList = {vStruct[], vArmour[]};
Field[5].IncludeBoundary = 1;

// F6-F8: FIRST WALL AND ARMOUR HELD FINE FROM THE HEATED FRONT TO CHANNEL 1.
Field[6] = Distance;
Field[6].SurfacesList = {sHeated[]};
Field[6].Sampling = 60;

Field[7] = Threshold;
Field[7].InField = 6;
Field[7].SizeMin = hFront;
Field[7].SizeMax = hSolid;
Field[7].DistMin = tArm + tFW;
Field[7].DistMax = tArm + tFW + dSolid;

Field[8] = Restrict;
Field[8].InField = 7;
Field[8].VolumesList = {vStruct[], vArmour[]};
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
Printf("BLANKET SLICE: TAPERED PLATE, 4 NESTED U-CHANNELS, UNIFORM OFFSETS");
Printf("back / front width, depth  = %.4e / %.4e / %.4e m", Wb, Wf, H);
Printf("plate front fillet R0      = %.4e m", R0);
Printf("slice thickness Tz         = %.4e m (z walls %.4e each side of channel 1)", Tz, tz);
Printf("armour / first wall        = %.4e / %.4e m", tArm, tFW);
Printf("gap between channels       = %.4e m; core cut-out %g at offset %.4e m (radius %.4e)",
       g, coreCut, dCore, R0 - dCore);
Printf("volumes fluid/struct/arm   = %g / %g / %g", #vFluid[], #vStruct[], #vArmour[]);
Printf("faces: inlets %g, outlets %g, heated %g, adiabatic %g, interface %g",
       #sInlet[], #sOutlet[], #sHeated[], #sAdiabatic[], #sInterface[]);
Printf("sizes: wall %.3e, core %.3e, front %.3e, solid %.3e m", hWall, hCore, hFront, hSolid);
Printf("==============================================================");
