// Gmsh project created on Tue Jun  4 10:03:57 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {-1, 1, 0, 1.0};
//+
Point(2) = {-0, 1, 0, 1.0};
//+
Point(3) = {-1, -0, 0, 1.0};
//+
Point(4) = {0, -0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 4};
//+
Line(3) = {3, 1};
//+
Line(4) = {4, 2};
//+
Curve Loop(1) = {3, 1, -4, -2};
//+
Plane Surface(1) = {1};
