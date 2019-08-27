// Gmsh project created on Tue Aug 20 09:50:45 2019
SetFactory("OpenCASCADE");

//+
Point(1) = {0, 1, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 0, 0, 1.0};
//+
Point(4) = {0.5, 0, 0, 1.0};
//+
Point(5) = {0, 0.5, 0, 1.0};
//+
Point(6) = {0.5, 0.5, 0, 1.0};
//+
Circle(1) = {1, 3, 2};
//+
Line(2) = {1, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {6, 4};
//+
Line(5) = {4, 2};
//+
Curve Loop(1) = {1, -5, -4, -3, -2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("dir_100") = {2};
//+
Physical Curve("dir_0") = {4};
//+
Physical Surface("air") = {1};
