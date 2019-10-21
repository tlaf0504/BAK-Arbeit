// Gmsh project created on Fri Oct 11 20:02:25 2019
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0.5, 0, 0.1};
//+
Point(2) = {0, -0.5, 0, 0.1};
//+
Point(3) = {1, -0.5, 0, 0.1};
//+
Point(4) = {1, 0.5, 0, 0.1};
//+
Circle(2) = {0.25, 0, 0, 0.005, 0, 2*Pi};
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 1};
//+
Curve Loop(1) = {6, 3, 4, 5};
//+
Curve Loop(2) = {2};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("dir0") = {6, 3, 5, 4};
//+
Physical Surface("air") = {1};
//+
Curve Loop(3) = {2};
//+
Plane Surface(2) = {3};
//+
Physical Surface("coil") = {2};
