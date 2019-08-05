// Gmsh project created on Fri Jul 19 17:47:42 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 10.0};
//+
Point(2) = {1, 0, 0, 10.0};
//+
Point(3) = {1, 0.5, 0, 10.0};
//+
Point(4) = {1, 1, 0, 10.0};
//+
Point(5) = {0, 0.5, 0, 10.0};
//+
Point(6) = {0, 1, 0, 10.0};
//+
Line(1) = {6, 5};
//+
Line(2) = {5, 1};
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 6};
//+
Line(7) = {5, 3};
//+
Curve Loop(1) = {6, 1, 7, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, -4, -3, -2};
//+
Plane Surface(2) = {2};
//+
Physical Surface("material2") = {2};
//+
Physical Surface("material1") = {1};
//+
Physical Curve("dir100") = {6};
//+
Physical Curve("dir0") = {3};
