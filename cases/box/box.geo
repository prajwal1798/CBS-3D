// Unit cube, free tetrahedra, three physical boundary groups and one fluid volume.
SetFactory("Built-in");

Point(1) = {0, 0, 0, 0.2};
Point(2) = {1, 0, 0, 0.2};
Point(3) = {1, 1, 0, 0.2};
Point(4) = {0, 1, 0, 0.2};
Point(5) = {0, 0, 1, 0.2};
Point(6) = {1, 0, 1, 0.2};
Point(7) = {1, 1, 1, 0.2};
Point(8) = {0, 1, 1, 0.2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Curve Loop(3) = {1, 10, -5, -9};
Plane Surface(3) = {3};

Curve Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4};

Curve Loop(5) = {3, 12, -7, -11};
Plane Surface(5) = {5};

Curve Loop(6) = {4, 9, -8, -12};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Names are convenience labels. The numeric tag is the important part.
Physical Surface("inlet_511", 511) = {6};
Physical Surface("outlet_520", 520) = {4};
Physical Surface("walls_530", 530) = {1, 2, 3, 5};
Physical Volume("fluid", 1) = {1};

Mesh.ElementOrder = 1;
Mesh.RecombineAll = 0;
Mesh.Recombine3DAll = 0;
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 4;
Mesh.Optimize = 1;
Mesh.Smoothing = 3;
Mesh.MshFileVersion = 2.2;
