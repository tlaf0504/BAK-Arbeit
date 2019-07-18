// Gmsh project created on Thu Jul 18 20:13:28 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {1, 0, 0, 0.1};
//+
Point(3) = {1, 1, 0, 0.1};
//+
Point(4) = {0, 1, 0, 0.1};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {3, 2};
//+
Line(4) = {3, 4};
//+
Circle(5) = {0.5, 0.5, 0, 0.15, 0, 2 * Pi};
//+
Curve Loop(1) = {4, 1, 2, -3};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};
//+
Physical Surface("air") = {1};
//+
Curve Loop(3) = {5};
//+
Plane Surface(2) = {3};
//+
Physical Surface("material") = {2};
//+
Physical Curve("dir100") = {4};
//+
Physical Curve("dir0") = {2};
