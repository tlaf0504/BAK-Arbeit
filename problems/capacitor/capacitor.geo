// Gmsh project created on Fri Mar 08 13:27:36 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.2};
//+
Point(2) = {1, 0, 0, .2};
//+
Point(3) = {1, 1, 0, .2};
//+
Point(4) = {0, 1, 0, .2};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Point(5) = {0.3, 0.52, 0, .2};
//+
Point(6) = {0.7, 0.52, 0, .2};
//+
Point(7) = {0.7, 0.53, 0, .2};
//+
Point(8) = {0.3, 0.53, 0, .2};
//+
Point(9) = {0.3, 0.48, 0,.2};
//+
Point(10) = {0.7, 0.48, 0, .2};
//+
Point(11) = {0.7, 0.47, 0, .2};
//+
Point(12) = {0.3, 0.47, 0, .2};
//+
Point(13) = {0.6, 0.8, 0, 1.0};
//+
Point(14) = {0.4, 0.8, 0, 1.0};
//+
Line(5) = {8, 7};
//+
Line(6) = {7, 6};
//+
Line(7) = {6, 5};
//+
Line(8) = {5, 8};
//+
Line(9) = {9, 12};
//+
Line(10) = {12, 11};
//+
Line(11) = {11, 10};
//+
Line(12) = {10, 9};
//+
Line(13) = {9, 5};
//+
Line(14) = {5, 6};
//+
Line(15) = {6, 10};
//+
Line(16) = {14, 13};


//+
Curve Loop(1) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, 11, 12, 9};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, -13, -12, -15};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 4, 1, 2};
//+
Curve Loop(5) = {13, 8, 5, 6, 15, -11, -10, -9};
//+
Plane Surface(4) = {4, 5};
//+
Physical Surface("air") = {4};
//+
Physical Surface("material") = {3};
//+
Physical Curve("dirichlet100") = {5, 8, 7, 6};
//+
Physical Curve("dirichlet0") = {12, 9, 10, 11};
//+
Physical Curve("farbound") = {4, 3, 2, 1};



//+
Physical Curve("neumann_1") = {16};
