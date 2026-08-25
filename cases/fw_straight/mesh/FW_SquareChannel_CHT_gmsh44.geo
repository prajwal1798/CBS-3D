//==============================================================================
// FW_SquareChannel_CHT_gmsh44.geo
//
// OpenCASCADE-free geometric equivalent of FW_SquareChannel_CHT.geo for the
// Sunbird gmsh/4.4.0 module, which is built without OpenCASCADE support.
//
// IMPORTANT NUMERICAL POINT:
//   This file does NOT use Extrude/Layers or a structured prism construction.
//   The geometry is assembled explicitly from points/lines/surfaces/volumes and
//   is still tetrahedralized by the genuine 3-D HXT algorithm (Algorithm3D=10).
//
// Geometry and physical topology are identical to the accepted first-wall case:
//   L                 = 200 mm
//   outer EUROFER     = 21 mm x 21 mm
//   helium channel    = 15 mm x 10 mm
//   fluid z           = 8 mm .. 18 mm
//   BC511             = fluid inlet only
//   BC520             = fluid outlet only
//   BC530             = all external solid surfaces except heated top
//   BC532             = external top surface
//   volume 10         = fluid
//   volume 20         = solid
//   internal interface= four shared surfaces, deliberately non-physical
//==============================================================================

//------------------------------------------------------------------------------
// 1. Geometry / mesh parameters
//------------------------------------------------------------------------------
L = 0.200;

Wout = 0.021;
Hout = 0.021;
Wfluid = 0.015;
Hfluid = 0.010;

tSide = 0.003;
tTop = 0.003;
tBottom = Hout - Hfluid - tTop;

yOuterMin = -Wout / 2.0;
yOuterMax =  Wout / 2.0;
zOuterMin = 0.0;
zOuterMax = Hout;

yFluidMin = -Wfluid / 2.0;
yFluidMax =  Wfluid / 2.0;
zFluidMin = tBottom;
zFluidMax = tBottom + Hfluid;

hInterface = 0.00015;
hFluidCore = 0.00060;
hBulk = 0.00100;
dFine = 0.00020;
dTransition = 0.00150;

// Characteristic lengths attached to points are intentionally neutral because
// explicit background fields below are the sizing authority.
lc = hBulk;

//------------------------------------------------------------------------------
// 2. Points: outer rectangle at x=0 and x=L
//------------------------------------------------------------------------------
Point(1) = {0.0, yOuterMin, zOuterMin, lc};
Point(2) = {0.0, yOuterMax, zOuterMin, lc};
Point(3) = {0.0, yOuterMax, zOuterMax, lc};
Point(4) = {0.0, yOuterMin, zOuterMax, lc};

Point(5) = {L, yOuterMin, zOuterMin, lc};
Point(6) = {L, yOuterMax, zOuterMin, lc};
Point(7) = {L, yOuterMax, zOuterMax, lc};
Point(8) = {L, yOuterMin, zOuterMax, lc};

// Inner helium rectangle at x=0 and x=L.
Point(9)  = {0.0, yFluidMin, zFluidMin, lc};
Point(10) = {0.0, yFluidMax, zFluidMin, lc};
Point(11) = {0.0, yFluidMax, zFluidMax, lc};
Point(12) = {0.0, yFluidMin, zFluidMax, lc};

Point(13) = {L, yFluidMin, zFluidMin, lc};
Point(14) = {L, yFluidMax, zFluidMin, lc};
Point(15) = {L, yFluidMax, zFluidMax, lc};
Point(16) = {L, yFluidMin, zFluidMax, lc};

//------------------------------------------------------------------------------
// 3. Curves
//------------------------------------------------------------------------------
// Outer x=0 perimeter.
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Outer x=L perimeter.
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Outer streamwise edges.
Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Inner x=0 perimeter.
Line(13) = {9, 10};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12, 9};

// Inner x=L perimeter.
Line(17) = {13, 14};
Line(18) = {14, 15};
Line(19) = {15, 16};
Line(20) = {16, 13};

// Inner streamwise edges.
Line(21) = {9, 13};
Line(22) = {10, 14};
Line(23) = {11, 15};
Line(24) = {12, 16};

//------------------------------------------------------------------------------
// 4. End-plane loops and surfaces
//------------------------------------------------------------------------------
// x=0 annular solid surface: outward normal -x.  The inner loop is reversed
// relative to the outer loop because it is a hole.
Line Loop(101) = {-4, -3, -2, -1};
Line Loop(102) = {13, 14, 15, 16};
Plane Surface(1001) = {101, 102};

// Fluid inlet at x=0: outward normal -x.
Line Loop(103) = {-16, -15, -14, -13};
Plane Surface(1101) = {103};

// x=L annular solid surface: outward normal +x.
Line Loop(104) = {5, 6, 7, 8};
Line Loop(105) = {-20, -19, -18, -17};
Plane Surface(1002) = {104, 105};

// Fluid outlet at x=L: outward normal +x.
Line Loop(106) = {17, 18, 19, 20};
Plane Surface(1102) = {106};

//------------------------------------------------------------------------------
// 5. Outer solid surfaces
//------------------------------------------------------------------------------
// Bottom z=zOuterMin, outward -z.
Line Loop(201) = {1, 10, -5, -9};
Plane Surface(1201) = {201};

// Right y=yOuterMax, outward +y.
Line Loop(202) = {2, 11, -6, -10};
Plane Surface(1202) = {202};

// Top z=zOuterMax, outward +z. This is the heat-flux surface.
Line Loop(203) = {3, 12, -7, -11};
Plane Surface(1203) = {203};

// Left y=yOuterMin, outward -y.
Line Loop(204) = {4, 9, -8, -12};
Plane Surface(1204) = {204};

//------------------------------------------------------------------------------
// 6. Four shared fluid-solid interface surfaces
//    Orientations below are outward from the FLUID volume.
//------------------------------------------------------------------------------
// Bottom interface, fluid outward -z.
Line Loop(301) = {13, 22, -17, -21};
Plane Surface(1301) = {301};

// Right interface, fluid outward +y.
Line Loop(302) = {14, 23, -18, -22};
Plane Surface(1302) = {302};

// Top interface, fluid outward +z.
Line Loop(303) = {15, 24, -19, -23};
Plane Surface(1303) = {303};

// Left interface, fluid outward -y.
Line Loop(304) = {16, 21, -20, -24};
Plane Surface(1304) = {304};

//------------------------------------------------------------------------------
// 7. Volumes
//------------------------------------------------------------------------------
// Fluid shell is directly closed by its inlet/outlet and the four interface
// surfaces.  All orientations are outward from fluid.
Surface Loop(2001) = {1101, 1102, 1301, 1302, 1303, 1304};
Volume(10) = {2001};

// Solid ring shell uses the two annular end surfaces, four external surfaces,
// and the REVERSED interface surfaces.  Reusing the exact same interface
// surfaces gives conformal fluid-solid topology without any OCC Boolean.
Surface Loop(2002) = {
    1001, 1002,
    1201, 1202, 1203, 1204,
    -1301, -1302, -1303, -1304
};
Volume(20) = {2002};

//------------------------------------------------------------------------------
// 8. CBS physical groups
//------------------------------------------------------------------------------
Physical Volume("fluid", 10) = {10};
Physical Volume("solid", 20) = {20};

Physical Surface("inlet", 511) = {1101};
Physical Surface("outlet", 520) = {1102};
Physical Surface("heat_flux", 532) = {1203};

// All remaining EXTERNAL solid surfaces, including both annular end faces.
Physical Surface("adiabatic_wall", 530) = {
    1001, 1002,
    1201, 1202, 1204
};

// Deliberately no Physical Surface for 1301..1304.  CBS3D reconstructs the CHT
// interface from shared topology + material adjacency.

//------------------------------------------------------------------------------
// 9. Mesh fields
//------------------------------------------------------------------------------
Field[1] = Distance;
Field[1].FacesList = {1301, 1302, 1303, 1304};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = hInterface;
Field[2].LcMax = hBulk;
Field[2].DistMin = dFine;
Field[2].DistMax = dTransition;
Field[2].Sigmoid = 1;

Field[3] = Box;
Field[3].VIn = hFluidCore;
Field[3].VOut = hBulk;
Field[3].XMin = 0.0;
Field[3].XMax = L;
Field[3].YMin = yFluidMin;
Field[3].YMax = yFluidMax;
Field[3].ZMin = zFluidMin;
Field[3].ZMax = zFluidMax;
Field[3].Thickness = 0.00050;

Field[4] = Min;
Field[4].FieldsList = {2, 3};
Background Field = 4;

// Gmsh-4.4 option names.
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthMin = hInterface;
Mesh.CharacteristicLengthMax = hBulk;

//------------------------------------------------------------------------------
// 10. Genuine 3-D linear tetrahedral mesh
//------------------------------------------------------------------------------
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 10;
Mesh.ElementOrder = 1;
Mesh.Smoothing = 10;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;
Mesh.SaveAll = 0;

//------------------------------------------------------------------------------
// 11. Audit output
//------------------------------------------------------------------------------
Printf("==============================================================");
Printf("FW200 OCC-FREE TRUE 3-D TET CHT MESH");
Printf("L                  = %.6e m", L);
Printf("Outer WxH          = %.6e x %.6e m", Wout, Hout);
Printf("Fluid WxH          = %.6e x %.6e m", Wfluid, Hfluid);
Printf("Fluid z range      = %.6e .. %.6e m", zFluidMin, zFluidMax);
Printf("Interface h target = %.6e m", hInterface);
Printf("Fluid-core target  = %.6e m", hFluidCore);
Printf("Solid bulk target  = %.6e m", hBulk);
Printf("Interface surfaces = 4");
Printf("BC511 inlet faces  = 1");
Printf("BC520 outlet faces = 1");
Printf("BC532 heat faces   = 1");
Printf("BC530 adiabatic    = 5");
Printf("==============================================================");
