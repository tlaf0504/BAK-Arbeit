// Gmsh project created on Tue Aug 20 08:23:51 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.01};
//+
Point(2) = {0.01, 0, 0, 0.01};
//+
Point(3) = {0.01, 0.001, 0, 0.01};
//+
Point(4) = {0, 0.001, 0, 0.01};
//+
Line(1) = {1, 2};
//+
Line(2) = {4, 1};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 2};
//+
Curve Loop(1) = {1, -4, -3, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("dir_10", 5) = {3};
//+
Physical Curve("dir_0", 6) = {1};
//+
Physical Surface("air", 7) = {1};
