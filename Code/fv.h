/*
 * fv.h
 *      Author: sunder
 */

#ifndef FV_H_
#define FV_H_

#include "triangulation.h"

/* CLS reconstruction on a face */

void clsq_construct_interior_face_data(int, const triangulation*, clsq_data*);
void clsq_construct_dirichlet_boundary_face_data(int, const triangulation*, clsq_data*);
void clsq_construct_neumann_boundary_face_data(int, const triangulation*, clsq_data*);

void assemble_gradient_term_interior_face(int, const triangulation*, const clsq_data*, double*);
void assemble_gradient_term_dirichlet_boundary_face(int, const triangulation*, const clsq_data*, double*);

/* Poisson equation solver. Solve the Poisson equation in the form -\nabla^2 u = f(x,y) */

enum poisson_solver_bcond { neumann, dirichlet };

typedef struct {

    triangulation* mesh;                // Mesh object associated with the solver
    double* phi;                        // Solution of the Poisson equation at cell centers of the mesh
    clsq_data* clsq;                    // Constrained least-squares data structure on each face
    double *phi_bnd;                    // Values of variable at boundary faces
    enum poisson_solver_bcond *bcond;   // Boundary condition array
    int *bmap;                          // Boundary condition map
    LIS_MATRIX A;                       // System matrix
    LIS_VECTOR f;                       // Forcing vector
    LIS_VECTOR u;                       // RHS vector
    LIS_SOLVER solver;                  // Iterative solver
} poisson_solver;


poisson_solver* poisson_solver_allocate(const char*);
void poisson_solver_free(poisson_solver*);
void poisson_solver_construct_clsq_data(poisson_solver*);
double poisson_solver_rhs(double,double);
double poisson_solver_exact_solution(double,double);
void poisson_solver_boundary_conditions(double, double, int, enum poisson_solver_bcond*,double*);
void poisson_solver_plot_vtk(poisson_solver*, const char*);
void poisson_solver_assemble_system(poisson_solver*);
void poisson_solver_run(poisson_solver*);

#endif /* FV_H_ */
