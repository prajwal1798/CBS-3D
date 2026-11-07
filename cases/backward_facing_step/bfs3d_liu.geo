// ============================================================================
// 3D TURBULENT BACKWARD-FACING STEP, Re = 3025 (LIU 2005, FIG. 7.6 AND SEC. 8.1)
//
// NON-DIMENSIONAL GEOMETRY, STEP HEIGHT L = 1:
//
//          y = 3  +------------------------------------------------------+
//    INLET        |                     TOP WALL                         |
//    (2L HIGH)    |                                                      | OUTLET
//          y = 1  +----------+                                           | p = 0
//                 0     x = 4 |  STEP FACE                               |
//          y = 0              +-------------------------------------------+
//                                        BOTTOM WALL                  x = 40
//
//   UPSTREAM CHANNEL  x = 0 .. 4,  y = 1 .. 3
//   DOWNSTREAM        x = 4 .. 40, y = 0 .. 3
//   SPAN              z = 0 .. 1   (THE 2D SECTION EXTRUDED BY UNIT LENGTH)
//
// BOUNDARY CODES (THE PHYSICAL-GROUP NUMBER IS THE CBS BOUNDARY CODE):
//   x = 0           INLET, MEASURED u(y), v = w = 0      -> 512 (PROFILE INLET)
//   x = 40          OUTLET, p = 0                        -> 520
//   ALL WALLS       NO-SLIP (TOP, UPSTREAM FLOOR, STEP FACE, DOWNSTREAM FLOOR)
//                                                        -> 530
//   z = 0, z = 1    SYMMETRY (2D-EQUIVALENT; LIU REPORTS IDENTICAL 2D AND 3D
//                   RESULTS, SO THE SIDE FACES DO NOT ACT AS WALLS)  -> 506
//
// MESH SIZE (FREE TETRAHEDRA, BACKGROUND FIELD; POINT SIZES ARE IGNORED):
//   lc_far     EVERYWHERE OUTSIDE THE RECIRCULATION BOX
//   lc_fine    INSIDE THE RECIRCULATION BOX x = x_box0 .. x_box1
//   lc_lip     AROUND THE STEP LIP (4, 1), WHERE THE SHEAR LAYER STARTS
//   ESTIMATE FOR THE DEFAULTS: ~ 0.3 M TETRAHEDRA, THE SAME ORDER AS LIU'S
//   297,054. GMSH PRINTS THE TRUE COUNT; MATCH IT WITH -clscale.
//   MESH STUDY: -clscale 1.0 (L1), 0.79 (L2, ~2x), 0.63 (L3, ~4x).
//
// MESH:   gmsh bfs3d_liu.geo -3 -format msh22 -o bfs3d_liu.msh
// ============================================================================
SetFactory("Built-in");

// ================================
// PARAMETERS
// ================================
L        = 1.0;     // STEP HEIGHT
x_step   = 4.0*L;   // INLET TO STEP
x_out    = 40.0*L;  // INLET TO OUTLET
h_in     = 2.0*L;   // UPSTREAM CHANNEL HEIGHT
span     = 1.0*L;   // EXTRUSION IN z

lc_far   = 0.25;
lc_fine  = 0.12;
lc_lip   = 0.05;
x_box0   = 3.5;
x_box1   = 16.0;

// ================================
// 2D SECTION AT z = 0
// ================================
Point(1) = {0.0,    L,        0.0, lc_far};   // INLET, FLOOR
Point(2) = {x_step, L,        0.0, lc_lip};   // STEP LIP
Point(3) = {x_step, 0.0,      0.0, lc_fine};  // STEP FOOT
Point(4) = {x_out,  0.0,      0.0, lc_far};   // OUTLET, FLOOR
Point(5) = {x_out,  L + h_in, 0.0, lc_far};   // OUTLET, ROOF
Point(6) = {0.0,    L + h_in, 0.0, lc_far};   // INLET, ROOF

Line(1) = {1, 2};   // UPSTREAM FLOOR
Line(2) = {2, 3};   // STEP FACE
Line(3) = {3, 4};   // DOWNSTREAM FLOOR
Line(4) = {4, 5};   // OUTLET
Line(5) = {5, 6};   // ROOF
Line(6) = {6, 1};   // INLET

Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// ================================
// EXTRUSION: out[0] = z = 1 FACE, out[1] = VOLUME,
// out[2..7] = SIDE FACES SWEPT FROM LINES 1..6 IN ORDER
// ================================
out[] = Extrude {0, 0, span} { Surface{1}; };

Physical Surface("inlet_profile_512", 512) = {out[7]};
Physical Surface("outlet_520", 520)        = {out[5]};
Physical Surface("walls_530", 530)         = {out[2], out[3], out[4], out[6]};
Physical Surface("symmetry_506", 506)      = {1, out[0]};
Physical Volume("fluid", 1)                = {out[1]};

// ================================
// SIZE FIELDS
// ================================
Field[1] = Box;
Field[1].VIn = lc_fine;
Field[1].VOut = lc_far;
Field[1].XMin = x_box0;  Field[1].XMax = x_box1;
Field[1].YMin = -1.0;    Field[1].YMax = L + h_in + 1.0;
Field[1].ZMin = -1.0;    Field[1].ZMax = span + 1.0;
Field[1].Thickness = 2.0;

Field[2] = Box;
Field[2].VIn = lc_lip;
Field[2].VOut = lc_far;
Field[2].XMin = x_step - 0.2;  Field[2].XMax = x_step + 0.6;
Field[2].YMin = L - 0.3;       Field[2].YMax = L + 0.3;
Field[2].ZMin = -1.0;          Field[2].ZMax = span + 1.0;
Field[2].Thickness = 0.6;

Field[3] = Min;
Field[3].FieldsList = {1, 2};
Background Field = 3;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6;          // 2D: FRONTAL-DELAUNAY
Mesh.Algorithm3D = 1;        // 3D: DELAUNAY
Mesh.Optimize = 1;
Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 2.2;
