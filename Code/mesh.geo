xmin = 0.0;
ymin = 0.0;
L = 1.;
np = 128.;

Point(1) = {xmin, ymin, 0.0, L/np};
Point(2) = {xmin+L, ymin, 0.0, L/np};
Point(3) = {xmin+L, ymin+L, 0.0, L/np};
Point(4) = {xmin, ymin+L, 0.0, L/np};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Recombine Surface {1};
Transfinite Surface {1};

Physical Curve("230") = {1,2,4}; // Left, right and bottom walls
Physical Curve("231") = {3};     // Top wall

// gmsh -2 -format su2 -save_all mesh.geo -o mesh.su2
