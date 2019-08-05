// Gmsh project created on Fri Aug  2 20:43:12 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {-20, -20, 0, 10};
//+
Point(2) = {20, -20, 0, 10};
//+
Point(3) = {20, 20, 0, 10};
//+
Point(4) = {-20, 20, 0, 10};
//+
Point(5) = {-0.5, -0.5, 0, 0.1};
//+
Point(6) = {0.5, -0.5, 0, 0.1};
//+
Point(7) = {0.5, 0.5, 0, 0.1};
//+
Point(8) = {-0.5, 0.5, 0, 0.1};
//+
Point(9) = {-0.5, -0.55, 0, 0.1};
//+
Point(10) = {0.5, -0.55, 0, 0.1};
//+
Point(11) = {0.5, 0.55, 0, 0.1};
//+
Point(12) = {-0.5, 0.55, 0, 0.1};
//+
Point(13) = {-0.5, 0, 0, 0.1};
//+
Point(14) = {0.5, 0, 0, 0.1};
//+
Line(1) = {12, 11};
//+
Line(2) = {11, 7};
//+
Line(3) = {8, 7};
//+
Line(4) = {12, 8};
//+
Line(5) = {5, 9};
//+
Line(6) = {9, 10};
//+
Line(7) = {10, 6};
//+
Line(8) = {6, 5};
//+
Line(9) = {8, 13};
//+
Line(10) = {13, 5};
//+
Line(11) = {7, 14};
//+
Line(12) = {14, 6};
//+
Line(13) = {14, 13};
//+
Line(14) = {4, 1};
//+
Line(15) = {1, 2};
//+
Line(16) = {2, 3};
//+
Line(17) = {3, 4};
//+
Physical Curve(1) = {17, 14, 15, 16};
//+
Physical Curve("farbound") = {14, 17, 16, 15};
//+
Physical Curve("dir100") = {1, 3, 4, 2};
//+
Physical Curve("dir0") = {8, 6, 5, 7};
//+
Curve Loop(1) = {9, -13, -11, -3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, -8, -12, 13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {16, 17, 14, 15};
//+
Curve Loop(4) = {11, 12, -7, -6, -5, -10, -9, -4, 1, 2};
//+
Plane Surface(3) = {3, 4};
//+
Physical Surface("material_1") = {1};
//+
Physical Surface("material_2") = {2};
//+
Physical Surface("air") = {3};
