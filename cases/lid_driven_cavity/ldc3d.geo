// ============================================================================
// 3D LID-DRIVEN CAVITY WITH HEATED SIDE WALLS - UNIT CUBE, FREE TETRAHEDRA
//
//   y = 1  : moving lid           -> CBS 500 (velocity = inlet_u v w in .par)
//   x = 0  : hot wall,  T = 1     -> CBS 501
//   x = 1  : cold wall, T = 0     -> CBS 502
//   y = 0, z = 0, z = 1 : walls   -> CBS 530 (no-slip, adiabatic)
//
// The physical-group NUMBER is the CBS boundary code, so no mapping is needed.
// Mesh size: change lc, or scale it on the command line with -clscale.
//   lc = 0.025  ->   ~0.29 M tets (old mesh)
//   lc = 0.0125 ->   ~2.3 M tets  (fine)
//   lc = 0.010  ->   ~4.5 M tets  (very fine; check memory first)
// ============================================================================
SetFactory("Built-in");

lc = 0.0125;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Point(5) = {0, 0, 1, lc};
Point(6) = {1, 0, 1, lc};
Point(7) = {1, 1, 1, lc};
Point(8) = {0, 1, 1, lc};

Line(1) = {1, 2};   Line(2) = {2, 3};   Line(3) = {3, 4};   Line(4) = {4, 1};
Line(5) = {5, 6};   Line(6) = {6, 7};   Line(7) = {7, 8};   Line(8) = {8, 5};
Line(9) = {1, 5};   Line(10) = {2, 6};  Line(11) = {3, 7};  Line(12) = {4, 8};

Curve Loop(1) = {1, 2, 3, 4};     Plane Surface(1) = {1};   // z = 0
Curve Loop(2) = {5, 6, 7, 8};     Plane Surface(2) = {2};   // z = 1
Curve Loop(3) = {1, 10, -5, -9};  Plane Surface(3) = {3};   // y = 0
Curve Loop(4) = {3, 12, -7, -11}; Plane Surface(4) = {4};   // y = 1 lid
Curve Loop(5) = {4, 9, -8, -12};  Plane Surface(5) = {5};   // x = 0 hot
Curve Loop(6) = {2, 11, -6, -10}; Plane Surface(6) = {6};   // x = 1 cold

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Physical Surface("walls_530", 530) = {1, 2, 3};
Physical Surface("lid_500", 500) = {4};
Physical Surface("hot_501", 501) = {5};
Physical Surface("cold_502", 502) = {6};
Physical Volume("fluid", 1) = {1};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.Algorithm = 6;          // 2D: Frontal-Delaunay
Mesh.Algorithm3D = 1;        // 3D: Delaunay (works in every Gmsh build)
Mesh.Optimize = 1;
Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 2.2;
