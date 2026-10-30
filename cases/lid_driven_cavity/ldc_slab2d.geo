// QUASI-2D LID-DRIVEN CAVITY: 1 x 1 x LZ SLAB, SYMMETRY (506) ON z = 0 AND z = LZ.
// WITH FRICTIONLESS FRONT/BACK FACES THE STEADY 3D SOLUTION IS THE 2D ONE,
// SO IT CAN BE COMPARED WITH GHIA ET AL. (1982) DIRECTLY.
SetFactory("Built-in");
lc = 0.0125;   // 80 x 80 IN-PLANE; 1.0/64 AND 1.0/128 FOR A MESH STUDY
Lz = 2*lc;
Point(1) = {0, 0, 0, lc}; Point(2) = {1, 0, 0, lc}; Point(3) = {1, 1, 0, lc}; Point(4) = {0, 1, 0, lc};
Point(5) = {0, 0, Lz, lc}; Point(6) = {1, 0, Lz, lc}; Point(7) = {1, 1, Lz, lc}; Point(8) = {0, 1, Lz, lc};
Line(1) = {1, 2};   Line(2) = {2, 3};   Line(3) = {3, 4};   Line(4) = {4, 1};
Line(5) = {5, 6};   Line(6) = {6, 7};   Line(7) = {7, 8};   Line(8) = {8, 5};
Line(9) = {1, 5};   Line(10) = {2, 6};  Line(11) = {3, 7};  Line(12) = {4, 8};
Curve Loop(1) = {1, 2, 3, 4};     Plane Surface(1) = {1};
Curve Loop(2) = {5, 6, 7, 8};     Plane Surface(2) = {2};
Curve Loop(3) = {1, 10, -5, -9};  Plane Surface(3) = {3};
Curve Loop(4) = {3, 12, -7, -11}; Plane Surface(4) = {4};
Curve Loop(5) = {4, 9, -8, -12};  Plane Surface(5) = {5};
Curve Loop(6) = {2, 11, -6, -10}; Plane Surface(6) = {6};
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Physical Surface("sym_506", 506) = {1, 2};
Physical Surface("walls_530", 530) = {3, 5, 6};
Physical Surface("lid_500", 500) = {4};
Physical Volume("fluid", 1) = {1};
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
Mesh.MshFileVersion = 2.2;
