lc = 0.01; // Define mesh characteristic length
R = 0.05;  // Radius of the cylinder
xc = 0.15+R;
yc = 0.15+R;

// Bounding box dimensions
xmin = 0.;
xmax = 2.2;
ymin = 0.;
ymax = 0.41;

// Define points (x, y, z, lc)
Point(1) = {xc, yc, 0, lc}; // Center
Point(2) = {xc-R, yc, 0, lc};
Point(3) = {xc, yc+R, 0, lc};
Point(4) = {xc+R, yc, 0, lc};
Point(5) = {xc, yc-R, 0, lc};

// Create circle arcs
Circle(1) = {2, 1, 3}; // StartPoint, CenterPoint, EndPoint
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Define Curve Loop and Surface
Curve Loop(1) = {1, 2, 3, 4};

Point(6) = {xmin,ymin,0,lc};
Point(7) = {xmax,ymin,0,lc};
Point(8) = {xmax,ymax,0,lc};
Point(9) = {xmin,ymax,0,lc};

Line(5) = {6,7};
Line(6) = {7,8};
Line(7) = {8,9};
Line(8) = {9,6};


Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {2, 1};

Delete { Point {1}; }



Physical Curve("9") = {8};       // Inlet
Physical Curve("23") = {1,2,3,4,5,7};  // Walls
Physical Curve("15") = {6};    // Outlet

//gmsh -2 -format su2 -save_all mesh.geo -o mesh.su2





