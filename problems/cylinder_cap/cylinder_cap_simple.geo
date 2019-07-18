// Gmsh project created on Mon Jun 17 15:59:56 2019
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 1, 0, Pi/2};
//+
Circle(2) = {0, 0, 0, 0.1, 0, Pi/2};
//+
Line(3) = {1, 3};
//+
Line(4) = {4, 2};
//+
Curve Loop(1) = {3, 2, 4, -1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("dir100") = {2};
//+
Physical Curve("dir0") = {1};
//+
Physical Surface("air") = {1};
//+
Physical Curve("neum") = {3};
