// Gmsh project created on Mon Jun 17 15:59:56 2019
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 1, 0, Pi/2};
//+
Circle(2) = {0, 0, 0, 0.1, 0, Pi/2};
//+
Circle(4) = {0, 0, 0, 0.5, 0, Pi/2};
//+
Line(5) = {1, 5};
//+
Line(6) = {5, 3};
//+
Line(7) = {4, 6};
//+
Line(8) = {6, 2};
//+
Curve Loop(1) = {6, 2, 7, -4};
//+
Plane Surface(1) = {1};
//+
Physical Surface("inner_material") = {1};
//+
Curve Loop(2) = {1, -8, -4, -5};
//+
Plane Surface(2) = {2};
//+
Physical Surface("outer_material") = {2};
//+
Physical Curve("dir100") = {1};
//+
Physical Curve("dir0") = {2};
