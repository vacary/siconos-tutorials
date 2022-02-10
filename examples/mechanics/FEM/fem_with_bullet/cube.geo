Point(1) = {0, 0., 0, 1.0};
Point(2) = {0, 0., 1, 1.0};
Point(3) = {0, 1., 0, 1.0};
Point(4) = {1, 0., 0, 1.0};
Point(5) = {1, 1., 0, 1.0};
Point(6) = {1, 0., 1, 1.0};
Point(7) = {0, 1., 1, 1.0};
Point(8) = {1, 1., 1, 1.0};



Line(1) = {1, 2};
Line(2) = {6, 2};
Line(3) = {6, 4};
Line(4) = {4, 1};
Line(5) = {1, 3};
Line(6) = {3, 7};
Line(7) = {7, 2};
Line(8) = {4, 5};
Line(9) = {5, 8};
Line(10) = {8, 6};
Line(11) = {7, 8};
Line(12) = {3, 5};
Line Loop(13) = {6, 11, -9, -12};
Plane Surface(14) = {13};
Line Loop(15) = {10, 3, 8, 9};
Plane Surface(16) = {15};
Line Loop(17) = {2, -1, -4, -3};
Plane Surface(18) = {17};
Line Loop(19) = {7, -1, 5, 6};
Plane Surface(20) = {19};
Line Loop(21) = {11, 10, 2, -7};
Plane Surface(22) = {21};
Line Loop(23) = {12, -8, 4, 5};
Plane Surface(24) = {23};
Surface Loop(25) = {22, 14, 20, 18, 24, 16};
Volume(26) = {25};
