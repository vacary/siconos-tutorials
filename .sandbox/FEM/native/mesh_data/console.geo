Point(1) = {-0.0, 0., 0, 1.0};
Point(2) = {10.0, 0., 0.75, 1.0};
Point(3) = {10.0, 1.0, 0.75, 1.0};
Point(4) = {10.0, 1.0, 1.0, 1.0};
Point(5) = {10.0, 0.0, 1.0, 1.0};
Line(1) = {5, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 5};
Point(6) = {0.0, 1.0, 0.0, 1.0};
Point(7) = {0.0, 1.0, 1.0, 1.0};
Point(8) = {0.0, 0.0, 1.0, 1.0};
Point(9) = {0.0, 0.0, 0.0, 1.0};
Line(5) = {7, 1};
Delete {
  Line{5};
}
Line(5) = {6, 1};
Line(6) = {6, 7};
Line(7) = {8, 8};
Line(8) = {1, 7};
Delete {
  Line{8};
}
Line(8) = {8, 1};
Line(9) = {7, 8};
Line(10) = {1, 2};
Line(11) = {3, 6};
Line(12) = {5, 8};
Line(13) = {4, 7};
Line Loop(14) = {5, 10, -3, 11};
Plane Surface(15) = {14};
Line Loop(16) = {8, 10, 4, 12};
Plane Surface(17) = {16};
Line Loop(18) = {6, 9, 8, -5};
Plane Surface(19) = {18};
Line Loop(20) = {11, 6, -13, 2};
Plane Surface(21) = {20};
Line Loop(22) = {1, 13, 9, -12};
Plane Surface(23) = {22};
Line Loop(24) = {4, 1, 2, 3};
Plane Surface(25) = {24};
Surface Loop(26) = {15, 19, 21, 23, 25, 17};
Volume(27) = {26};
