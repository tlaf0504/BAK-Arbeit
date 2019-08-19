// Gmsh project created on Sat Aug 17 12:31:48 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {-2.5, 0, 0, 1.0};
//+
Point(2) = {2.5, 0, 0, 1.0};
//+
Point(3) = {2.5, 5, 0, 1.0};
//+
Point(4) = {-2.5, 5, 0, 1.0};
//+
Point(5) = {0,0,0,0.01};
//+
Disk(1) = {-0.05, 0.1, 0, 0.02, 0.02};
//+
Disk(2) = {0.05, 0.2, 0, 0.02, 0.02};
//+
Line(3) = {1, 5};
//+
Line(4) = {5, 2};
//+
Line(5) = {2, 3};
//+
Line(6) = {3, 4};
//+
Line(7) = {4, 1};
//+
Curve Loop(3) = {6, 7, 3, 4, 5};
//+
Curve Loop(4) = {2};
//+
Curve Loop(5) = {1};
//+
Plane Surface(3) = {3, 4, 5};
//+
Physical Curve("dir_electrode_1") = {1};
//+
Physical Curve("dir_electrode_2") = {2};
//+
Physical Curve("farbound") = {7, 6, 5, 4, 3};
//+
Physical Surface("air") = {3};
