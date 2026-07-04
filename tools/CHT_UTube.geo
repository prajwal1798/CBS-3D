// ============================================================================
//  CHT_UTube.geo  --  U-tube coolant manifold in a RAFM-steel block
//  Conjugate Heat Transfer benchmark for CBS3D++_SI
//
//  Units: METRES (SI).  Geometry from SpaceClaim drawing:
//    Block 300 x 100 x 100 mm ; pipe bore D20 mm ; legs 260 mm ;
//    U-bend mean radius 25 mm (R15 inner / R35 outer).
//
//  Two conformal volumes share the pipe-wall interface (BooleanFragments),
//  so the conjugate flux continuity the solver needs is exact (shared nodes).
//
//  Flow  : inlet & outlet on the x = 0 face (mass-flow inlet / pressure outlet)
//  Thermal: heat flux on the far end face x = L (opposite the ports)
//
//  Build a mesh:   gmsh CHT_UTube.geo -3 -format msh2 -o CHT_UTube.msh
//  then convert .msh -> solver .plt/.bco with your GMSH->CBS converter,
//  using the PHYSICAL IDS documented at the bottom of this file.
// ============================================================================

SetFactory("OpenCASCADE");

// ---- parameters (metres) ---------------------------------------------------
L   = 0.300;     // block length  (x)
H   = 0.100;     // block height  (y)
W   = 0.100;     // block width   (z)

r    = 0.010;    // pipe radius (D20 mm)
Lleg = 0.260;    // straight-leg length to bend centre (x)
Rb   = 0.025;    // U-bend mean centreline radius (R15 inner, R35 outer)

yc   = 0.050;    // pipe centre height (mid-block)
zc1  = 0.025;    // INLET  leg centreline z
zc2  = 0.075;    // OUTLET leg centreline z   (zc2-zc1 = 2*Rb = 0.050)
zmid = 0.050;    // bend axis z

// ---- mesh sizes (metres) ---------------------------------------------------
lc_fine   = 0.0020;   // near pipe wall / interface, inlet, heat-flux face
lc_coarse = 0.0120;   // bulk solid
ref_dmin  = 0.0030;   // full refinement within this distance of refined faces
ref_dmax  = 0.0250;   // coarsen out to here

eps = 1e-5;           // bounding-box tolerance

// ============================================================================
//  1. FLUID U-TUBE
// ============================================================================
// straight legs (base centre, axis vector, radius)
Cylinder(1) = {0, yc, zc1, Lleg, 0, 0, r};   // inlet  leg
Cylinder(2) = {0, yc, zc2, Lleg, 0, 0, r};   // outlet leg

// U-bend: revolve the inlet-leg end disk 180 deg about the y-axis at the bend
// centre. Disk(101) is created normal to +z, so first rotate it to face +x.
Disk(101) = {Lleg, yc, zc1, r};
Rotate { {0,1,0}, {Lleg, yc, zc1}, Pi/2 } { Surface{101}; }
// -Pi makes the bend bulge toward +x (away from the inlet face).
// >>> If a first mesh shows the bend curving back over the legs, change -Pi to +Pi. <<<
bend[] = Extrude { {0,1,0}, {Lleg, yc, zmid}, -Pi } { Surface{101}; };

// fuse the three pieces into a single fluid solid
fluidall[] = BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Volume{bend[1]}; Delete; };

// ============================================================================
//  2. SOLID BLOCK  +  conformal split
// ============================================================================
Box(200) = {0, 0, 0, L, H, W};

// Fragment block against fluid -> two conformal volumes sharing the interface.
all[] = BooleanFragments{ Volume{200}; Delete; }{ Volume{fluidall[0]}; Delete; };

// ---- identify fluid vs solid by bounding box -------------------------------
// The fluid is thin in y and z (the block fills the full 0.1 there), so we
// query a box that spans the full x-length but is tight in y/z. The block's
// bbox is NOT inside this (it spans all y,z) so only the fluid is returned.
// pad absorbs the small margin OCC adds to curved-surface bounding boxes.
pad = 0.005;
fluidVol[] = Volume In BoundingBox{
    -eps,        yc-r-pad,  zc1-r-pad,
    L+eps,       yc+r+pad,  zc2+r+pad };
// solid = everything else
solidVol[] = {};
For i In {0 : #all[]-1}
  isFluid = 0;
  For j In {0 : #fluidVol[]-1}
    If (all[i] == fluidVol[j])
      isFluid = 1;
    EndIf
  EndFor
  If (isFluid == 0)
    solidVol[] += all[i];
  EndIf
EndFor

// ============================================================================
//  3. BOUNDARY SURFACES (located by bounding box -> robust to renumbering)
// ============================================================================
// INLET  : fluid disk on x=0 around z=zc1
inlet[] = Surface In BoundingBox{
    -eps, yc-r-pad, zc1-r-pad,  +eps, yc+r+pad, zc1+r+pad };
// OUTLET : fluid disk on x=0 around z=zc2
outlet[] = Surface In BoundingBox{
    -eps, yc-r-pad, zc2-r-pad,  +eps, yc+r+pad, zc2+r+pad };

// HEAT-FLUX face : entire far end face x=L
heatflux[] = Surface In BoundingBox{
    L-eps, -eps, -eps,  L+eps, H+eps, W+eps };

// PIPE WALL / INTERFACE : all fluid-volume faces except inlet & outlet
fb[] = Abs(Boundary{ Volume{fluidVol[0]}; });
wall[] = {};
For i In {0 : #fb[]-1}
  keep = 1;
  For j In {0 : #inlet[]-1}
    If (fb[i] == inlet[j])
      keep = 0;
    EndIf
  EndFor
  For j In {0 : #outlet[]-1}
    If (fb[i] == outlet[j])
      keep = 0;
    EndIf
  EndFor
  If (keep == 1)
    wall[] += fb[i];
  EndIf
EndFor

// ADIABATIC : all remaining outer block faces (everything not flux/inlet/outlet)
// = boundary of the solid volume minus the interface(=wall) minus heat-flux face.
sb[] = Abs(Boundary{ Volume{solidVol[0]}; });
adia[] = {};
For i In {0 : #sb[]-1}
  keep = 1;
  For j In {0 : #wall[]-1}
    If (sb[i] == wall[j])
      keep = 0;
    EndIf
  EndFor
  For j In {0 : #heatflux[]-1}
    If (sb[i] == heatflux[j])
      keep = 0;
    EndIf
  EndFor
  If (keep == 1)
    adia[] += sb[i];
  EndIf
EndFor

// ============================================================================
//  4. MESH REFINEMENT  (distance from interface, inlet, heat-flux face)
// ============================================================================
Field[1] = Distance;
Field[1].SurfacesList = { wall[], inlet[], heatflux[] };   // gmsh >= 4.7
// (older gmsh: use  Field[1].FacesList = { wall[], inlet[], heatflux[] };)

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_fine;
Field[2].SizeMax  = lc_coarse;
Field[2].DistMin  = ref_dmin;
Field[2].DistMax  = ref_dmax;

Background Field = 2;

Mesh.CharacteristicLengthMin = lc_fine;
Mesh.CharacteristicLengthMax = lc_coarse;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 12;   // a few elements around the pipe circumference
Mesh.Algorithm3D = 1;              // 1=Delaunay (robust); try 10=HXT for speed
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

// ============================================================================
//  5. PHYSICAL GROUPS   (IDs = raw flags carried into the .msh)
// ============================================================================
//  Map these raw flags to solver bc_codes in the .bco file:
//    501 inlet  -> mass-flow / velocity inlet   (e.g. solver code 511)
//    502 outlet -> pressure outlet              (e.g. solver code 504)
//    530 wall   -> no-slip wall + conjugate interface
//    531 adiab  -> adiabatic (natural / zero-flux); no velocity BC on solid faces
//    532 flux   -> applied heat flux 0.5 MW/m^2 = 5e5 W/m^2
//  Volumes:
//    10 fluid   -> mat_elem = fluid  (water properties)
//    20 solid   -> mat_elem = solid  (RAFM steel properties)
// ----------------------------------------------------------------------------
Physical Surface("inlet",     501) = { inlet[] };
Physical Surface("outlet",    502) = { outlet[] };
Physical Surface("wall",      530) = { wall[] };
Physical Surface("adiabatic", 531) = { adia[] };
Physical Surface("heatflux",  532) = { heatflux[] };

Physical Volume("fluid", 10) = { fluidVol[] };
Physical Volume("solid", 20) = { solidVol[] };
