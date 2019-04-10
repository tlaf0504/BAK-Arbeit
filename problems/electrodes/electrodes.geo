// Gmsh project created on Wed Feb 20 09:54:45 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 6, 0, 1.0};
//+
Point(4) = {0, 6, 0, 1.0};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Circle(5) = {3, 3, 0, 0.3, 0, 2*Pi};
//+
Circle(6) = {7, 3, 0, 0.3, 0, 2*Pi};
//+
Physical Curve("dirichlet_0") = {2};
//+
Curve Loop(1) = {5};
//+
Plane Surface(1) = {1};
//+
Physical Surface("source1") = {1};
//+
Curve Loop(2) = {6};
//+
Plane Surface(2) = {2};
//+
Physical Surface("source2") = {2};
//+
Curve Loop(3) = {4, 1, 2, 3};
//+
Plane Surface(3) = {3,1,2};
//+
Physical Surface("air") = {3};