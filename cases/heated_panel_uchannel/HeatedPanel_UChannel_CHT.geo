//==============================================================================
// HeatedPanel_UChannel_CHT.geo
//
// GENERIC CONJUGATE-HEAT-TRANSFER UNIT CELL FOR CBS3D++_SI:
// A SOLID PANEL HEATED ON ONE FACE AND COOLED BY A SQUARE U-CHANNEL.
//
// The cell is deliberately generic. It represents the common building block
// of actively cooled heated walls (first-wall / blanket panels, cold plates,
// turbine-blade internal cooling passages): thin wall between the heated face
// and the coolant, two parallel passes joined by a 180-degree bend, and an
// unheated entrance length. No specific blanket design is copied.
//
//                    heated face (z = T), x >= xHeat : BC 532
//        +-----------------------------------------------------------+
//        |  unheated  |                                               |
//   y    |  =======================================================\  |
//   ^  ->| inlet leg (y = -p/2)                                     |  |
//   |    |                                                          |  | bend, mean
//   +->x |  outlet leg (y = +p/2)                                   |  | radius p/2
//      <-|  =======================================================/  |
//        +-----------------------------------------------------------+
//        x = 0 : inlet (BC 511) and outlet (BC 520) ports, rest adiabatic
//
//   cross-section (x = const):        z ^
//        +-----------------------+    T |  - - heated face - -
//        |       tFront          |      |
//        |   +---+       +---+   |      |   channel a x a, centre z = zc
//        |   | a |       | a |   |      |
//        |   +---+       +---+   |      |
//        |       tBack           |    0 +-----------------------> y
//        +-----------------------+         -W/2      0      +W/2
//
// Units: METRES. Flow enters the leg at y = -p/2 and leaves at y = +p/2.
//
// CBS physical IDs (read by tools/gmsh_to_cbs.py):
//    10 fluid volume              20 solid volume
//   511 inlet (mass flow + T)    520 pressure outlet
//   532 heated face (flux)       530 adiabatic external solid surfaces
//   The fluid-solid interface is conformal (shared nodes) and is NOT a
//   physical boundary: the solver couples it through the shared nodes.
//
// Mesh: genuine 3-D linear tetrahedra (no extrusion, no prisms). Fluid and
// solid are graded separately: fine at the interface, a fine fluid core,
// a refined heated wall, and a coarser solid bulk.
//
// Every parameter can be changed without editing this file, e.g.
//   gmsh HeatedPanel_UChannel_CHT.geo -3 -setnumber a 0.008 -setnumber hScale 2
// hScale multiplies every mesh size (2 = about 8x fewer elements).
//==============================================================================

SetFactory("OpenCASCADE");
Geometry.OCCBooleanPreserveNumbering = 1;

//==============================================================================
// 1. GEOMETRY PARAMETERS
//==============================================================================

DefineConstant[
    a      = {0.010, Name "Geometry/1 channel side a (m)"},
    p      = {0.025, Name "Geometry/2 leg centre spacing p (m)"},
    Lleg   = {0.250, Name "Geometry/3 straight leg length (m)"},
    tFront = {0.003, Name "Geometry/4 heated wall thickness (m)"},
    tBack  = {0.007, Name "Geometry/5 back wall thickness (m)"},
    mEnd   = {0.005, Name "Geometry/6 solid beyond the bend (m)"},
    xHeat  = {0.050, Name "Geometry/7 start of heated face (m), 0 = all"}
];

// DERIVED DIMENSIONS
Rb = p / 2.0;                   // MEAN BEND RADIUS (LEG CENTRE LINES MEET)
W  = 2.0 * p;                   // CELL WIDTH: ONE U-LOOP PER CELL
T  = tBack + a + tFront;        // PANEL THICKNESS
zc = tBack + a / 2.0;           // CHANNEL CENTRE HEIGHT
Lp = Lleg + Rb + a / 2.0 + mEnd;// PANEL LENGTH

yIn  = -p / 2.0;                // INLET LEG CENTRE
yOut =  p / 2.0;                // OUTLET LEG CENTRE

eps = 1.0e-6;                   // ENTITY-SELECTION TOLERANCE ONLY

If (p <= a)
    Error("p (%g) must exceed a (%g): the bend inner radius p/2 - a/2 must be positive", p, a);
    Abort;
EndIf
If (xHeat < 0 || xHeat >= Lleg)
    Error("xHeat (%g) must lie in [0, Lleg)", xHeat);
    Abort;
EndIf

//==============================================================================
// 2. FLUID: TWO STRAIGHT LEGS + 180-DEGREE BEND (SQUARE SECTION)
//==============================================================================

vLegIn = newv;
Box(vLegIn) = {0.0, yIn - a / 2.0, zc - a / 2.0, Lleg, a, a};

vLegOut = newv;
Box(vLegOut) = {0.0, yOut - a / 2.0, zc - a / 2.0, Lleg, a, a};

// SQUARE SECTION AT THE END OF THE INLET LEG: BUILT IN THE XY PLANE, TURNED
// INTO THE YZ PLANE, MOVED TO x = Lleg, THEN SWEPT 180 DEGREES ABOUT THE
// z-AXIS THROUGH (Lleg, 0). +Pi SWEEPS THROUGH x > Lleg (AWAY FROM THE PORTS).
sSection = news;
Rectangle(sSection) = {0.0, 0.0, 0.0, a, a};
Rotate { {0, 1, 0}, {0, 0, 0}, -Pi / 2.0 } { Surface{sSection}; }
Translate { Lleg, yIn - a / 2.0, zc - a / 2.0 } { Surface{sSection}; }
bend[] = Extrude { {0, 0, 1}, {Lleg, 0, 0}, Pi } { Surface{sSection}; };

vFluidRaw[] = BooleanUnion { Volume{vLegIn}; Delete; } { Volume{vLegOut, bend[1]}; Delete; };

//==============================================================================
// 3. SOLID PANEL (SPLIT AT x = xHeat SO THE HEATED FACE STARTS THERE)
//==============================================================================

vPanel[] = {};
If (xHeat > 0)
    vA = newv;
    Box(vA) = {0.0,   -W / 2.0, 0.0, xHeat,      W, T};
    vB = newv;
    Box(vB) = {xHeat, -W / 2.0, 0.0, Lp - xHeat, W, T};
    vPanel[] = {vA, vB};
Else
    vA = newv;
    Box(vA) = {0.0, -W / 2.0, 0.0, Lp, W, T};
    vPanel[] = {vA};
EndIf

// CONFORMAL SPLIT: SOLID AND FLUID SHARE THE INTERFACE NODES.
allVol[] = BooleanFragments { Volume{vPanel[]}; Delete; } { Volume{vFluidRaw[]}; Delete; };

//==============================================================================
// 4. FLUID / SOLID VOLUMES (BY BOUNDING BOX)
//==============================================================================

// THE CHANNEL REGION IS THINNER THAN THE PANEL IN z, SO ONLY FLUID PIECES FIT.
pad = 1.0e-3;
vFluid[] = Volume In BoundingBox {
    -eps,                yIn - a / 2.0 - pad,  zc - a / 2.0 - pad,
    Lleg + Rb + a / 2.0 + pad,
                         yOut + a / 2.0 + pad, zc + a / 2.0 + pad };

vSolid[] = allVol[];
vSolid[] -= vFluid[];

If (#vFluid[] < 1 || #vSolid[] < 1)
    Error("Volume identification failed: %g fluid, %g solid volumes", #vFluid[], #vSolid[]);
    Abort;
EndIf

Physical Volume("fluid", 10) = {vFluid[]};
Physical Volume("solid", 20) = {vSolid[]};

//==============================================================================
// 5. BOUNDARY SURFACES
//==============================================================================

sInlet[] = Surface In BoundingBox {
    -eps, yIn - a / 2.0 - eps, zc - a / 2.0 - eps,
     eps, yIn + a / 2.0 + eps, zc + a / 2.0 + eps };

sOutlet[] = Surface In BoundingBox {
    -eps, yOut - a / 2.0 - eps, zc - a / 2.0 - eps,
     eps, yOut + a / 2.0 + eps, zc + a / 2.0 + eps };

sHeated[] = Surface In BoundingBox {
    xHeat - eps, -W / 2.0 - eps, T - eps,
    Lp + eps,     W / 2.0 + eps, T + eps };

// EVERY EXTERNAL SURFACE THAT IS NOT A PORT OR THE HEATED FACE IS ADIABATIC.
sExternal[] = Abs(CombinedBoundary { Volume{allVol[]}; });
sAdiabatic[] = sExternal[];
sAdiabatic[] -= sInlet[];
sAdiabatic[] -= sOutlet[];
sAdiabatic[] -= sHeated[];

If (#sInlet[] != 1 || #sOutlet[] != 1 || #sHeated[] < 1)
    Error("Boundary identification failed: inlet %g, outlet %g, heated %g",
          #sInlet[], #sOutlet[], #sHeated[]);
    Abort;
EndIf

Physical Surface("inlet", 511)          = {sInlet[]};
Physical Surface("outlet", 520)         = {sOutlet[]};
Physical Surface("heat_flux", 532)      = {sHeated[]};
Physical Surface("adiabatic_wall", 530) = {sAdiabatic[]};

// FLUID-SOLID INTERFACE (MESH SIZING ONLY): FLUID BOUNDARY MINUS PORTS MINUS
// THE INTERNAL FLUID CUT AT x = xHeat.
sInterface[] = Abs(Boundary { Volume{vFluid[]}; });
sInterface[] -= sInlet[];
sInterface[] -= sOutlet[];
If (xHeat > 0)
    sCut[] = Surface In BoundingBox {
        xHeat - eps, -W / 2.0 - eps, -eps,
        xHeat + eps,  W / 2.0 + eps, T + eps };
    sInterface[] -= sCut[];
EndIf

//==============================================================================
// 6. MESH SIZES
//==============================================================================

DefineConstant[
    hScale = {1.0,     Name "Mesh/0 global size factor (2 = coarser)"},
    hWall  = {0.15e-3, Name "Mesh/1 size at fluid-solid interface (m)"},
    hCore  = {0.45e-3, Name "Mesh/2 size in fluid core (m)"},
    hFront = {0.50e-3, Name "Mesh/3 size in heated wall (m)"},
    hSolid = {1.20e-3, Name "Mesh/4 size in solid bulk (m)"},
    dWall  = {0.10e-3, Name "Mesh/5 keep hWall up to this wall distance (m)"},
    dCore  = {2.00e-3, Name "Mesh/6 fluid reaches hCore at this distance (m)"},
    dSolid = {4.00e-3, Name "Mesh/7 solid reaches hSolid at this distance (m)"}
];

// FLUID AND SOLID ARE GRADED SEPARATELY (Restrict FIELDS): THE FLUID GROWS
// SMOOTHLY FROM hWall TO A FINE CORE hCore ACROSS THE CHANNEL; THE SOLID GROWS
// FASTER TO hSolid, WHERE CONDUCTION NEEDS NO FINE CELLS.

// F1: DISTANCE TO THE FLUID-SOLID INTERFACE.
Field[1] = Distance;
Field[1].SurfacesList = {sInterface[]};
Field[1].Sampling = 80;

// F2-F3: FLUID GRADATION, ACTIVE INSIDE THE FLUID ONLY.
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = hScale * hWall;
Field[2].SizeMax = hScale * hCore;
Field[2].DistMin = dWall;
Field[2].DistMax = dCore;

Field[3] = Restrict;
Field[3].InField = 2;
Field[3].VolumesList = {vFluid[]};
Field[3].IncludeBoundary = 1;

// F4-F5: SOLID GRADATION, ACTIVE INSIDE THE SOLID ONLY.
Field[4] = Threshold;
Field[4].InField = 1;
Field[4].SizeMin = hScale * hWall;
Field[4].SizeMax = hScale * hSolid;
Field[4].DistMin = dWall;
Field[4].DistMax = dSolid;

Field[5] = Restrict;
Field[5].InField = 4;
Field[5].VolumesList = {vSolid[]};
Field[5].IncludeBoundary = 1;

// F6: HEATED WALL (BETWEEN THE CHANNEL TOP AND THE HEATED FACE): AT LEAST
// SIX CELLS THROUGH tFront FOR THE STEEPEST SOLID TEMPERATURE GRADIENT.
Field[6] = Box;
Field[6].VIn = hScale * hFront;
Field[6].VOut = hScale * hSolid;
Field[6].XMin = -eps;
Field[6].XMax = Lp + eps;
Field[6].YMin = -W / 2.0 - eps;
Field[6].YMax =  W / 2.0 + eps;
Field[6].ZMin = zc + a / 2.0;
Field[6].ZMax = T + eps;
Field[6].Thickness = tFront;

Field[7] = Min;
Field[7].FieldsList = {3, 5, 6};
Background Field = 7;

// ONLY THE BACKGROUND FIELD CONTROLS THE SIZE.
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

//==============================================================================
// 7. TRUE 3-D LINEAR TETRAHEDRAL MESH
//==============================================================================

DefineConstant[
    algo3D = {10, Name "Mesh/7 3-D algorithm (10 HXT, 1 Delaunay)"}
];

Mesh.Algorithm = 6;              // 2-D: FRONTAL-DELAUNAY
Mesh.Algorithm3D = algo3D;       // 3-D: HXT (PARALLEL), FALL BACK TO 1
Mesh.ElementOrder = 1;           // TET4 FOR CBS3D
Mesh.Smoothing = 10;
Mesh.Optimize = 1;               // GMSH TETRAHEDRAL OPTIMISER
Mesh.OptimizeNetgen = 1;         // NETGEN OPTIMISER (NEEDS A FULL GMSH BUILD)
Mesh.OptimizeThreshold = 0.35;   // IMPROVE EVERY TET BELOW THIS QUALITY
Mesh.MshFileVersion = 2.2;       // FORMAT READ BY tools/gmsh_to_cbs.py
Mesh.Binary = 0;
Mesh.SaveAll = 0;

//==============================================================================
// 8. AUDIT
//==============================================================================

Printf("==============================================================");
Printf("HEATED PANEL WITH U-CHANNEL: CHT UNIT CELL");
Printf("channel a x a           = %.4e m", a);
Printf("leg spacing p, bend Rb  = %.4e m, %.4e m (Rb/a = %.3f)", p, Rb, Rb / a);
Printf("panel L x W x T         = %.4e x %.4e x %.4e m", Lp, W, T);
Printf("heated wall / back wall = %.4e / %.4e m", tFront, tBack);
Printf("heated from x           = %.4e m", xHeat);
Printf("fluid volumes / solid   = %g / %g", #vFluid[], #vSolid[]);
Printf("faces: inlet %g, outlet %g, heated %g, adiabatic %g, interface %g",
       #sInlet[], #sOutlet[], #sHeated[], #sAdiabatic[], #sInterface[]);
Printf("sizes: wall %.3e, core %.3e, heated wall %.3e, solid %.3e (x hScale %g)",
       hWall, hCore, hFront, hSolid, hScale);
Printf("==============================================================");
