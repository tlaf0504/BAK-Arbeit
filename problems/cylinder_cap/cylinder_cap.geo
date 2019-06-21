// Gmsh project created on Mon Jun 17 15:59:56 2019
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 1, 0, Pi/2};
//+
Circle(2) = {0, 0, 0, 0.1, 0, Pi/2};
//+
Circle(3) = {0, 0, 0, 0.5, 0, Pi/2};
//+
Line(4) = {3, 5};
//+
Line(5) = {5, 1};
//+
Line(6) = {4, 6};
//+
Line(7) = {6, 2};
//+
Curve Loop(1) = {5, 1, -7, -3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, 3, -6, -2};
//+
Plane Surface(2) = {2};
//+
Physical Surface("material_1") = {2};
//+
Physical Surface("material_2") = {1};
//+
Physical Curve("dir0") = {1};
//+
Physical Curve("dir100") = {2};
//+
Physical Curve("neumann_1") = {3};
