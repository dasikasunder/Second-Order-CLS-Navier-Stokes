/*
 * triangulation.h
 *      Author: sunder
 */

#ifndef TRIANGULATION_H_
#define TRIANGULATION_H_

#include "includes.h"

/* Vertex data type */

typedef struct {
    double x[2]; // Coordinates of the vertex
    int nc;      // No. of cells surrounding the vertex
    int *c;      // List of cells surrounding the vertex
} Vrtx;

/* Cell data type */

typedef struct {
    int nv;       // No. of vertices in the cell
    int nf;       // No. of faces in the cell
    int v[4];     // List of vertices
    int f[4];     // List of faces
    int nn;       // No. of (face) neighbours of the cell
    int n[4];     // List of face neighbours (-1 if neighbour is not present)
    double sn[4]; // Sign of face normal
    double x[2];  // Cell centroid
    double vol;   // Volume of the cell (area in 2D)
} Cell;

/* Face data type */

typedef struct {
    int v[2];         // Vertices of the face
    int c[2];         // Cells adjacent to the face (for boundary faces c[1] = -1)
    int btag;         // Boundary tag associated with the face (0 for interior faces)
    double area;      // Area of the face (length in 2D)
    double n[2];      // Face normal
    double e[2];      // Unit vector connecting adjacent cell centers
    double r;         // Distance between adjacent cell centers
    double x[2];      // Center of the face
    int nvnb;         // No. of vertex neighbours (c[0] and c[1] are not counted in this list)
    int* vnb;         // List of vertex neighbours
    bool at_boundary; // Whether the face is at boundary
} Face;

/* Constrained least-squares data type */

typedef struct {
    double **sol; //
} clsq_data;


/* Triangulation object */

typedef struct  {
    int n_dim;          // No. of dimensions in the mesh (Must be two)
    int n_cell;         // No. of cells in the mesh
    int n_vrtx;         // No. of vertices in the mesh
    int n_face;         // No. of faces in the mesh
    int n_tria;         // No. of triangles in the mesh
    int n_quad;         // No. of quadrilaterals in the mesh
    int n_bnd_tag;      // No. of boundary tags in the mesh
    int n_int_face;     // No. of interior face
    int n_bnd_face;     // No. of faces on the boundary
    Vrtx* vrtx;         // Array of vertices in the mesh
    Cell* cell;         // Array of cells in the mesh
    Face* face;         // Array of faces in the mesh

} triangulation;

double triangle_area(double, double, double, double, double, double);
triangulation* triangulation_allocate(const char*);
void triangulation_plot_vtk(const triangulation*, const char*);
void triangulation_verify(const triangulation*);
void triangulation_free(triangulation*);

/* Cell-centered least squares data type */

typedef struct {
	double *ai;
    double **ak;
} cclsq_data;

void cclsq_data_initialize(cclsq_data*, const triangulation*);
void cclsq_data_free(cclsq_data*, int);

#endif /* TRIANGULATION_H_ */
