/*
 * fv.h
 *      Author: sunder
 */

#ifndef FV_H_
#define FV_H_

#include "triangulation.h"

/* Data type for storing reconstruction */

typedef struct {
    double **sol; //
} recon_data;

/* Constrained least-squares reconstruction on a face */

void clsq_construct_interior_face_data(int, const triangulation*, recon_data*);
void clsq_construct_dirichlet_boundary_face_data(int, const triangulation*, recon_data*);
void clsq_construct_neumann_boundary_face_data(int, const triangulation*, recon_data*);

void assemble_gradient_term_interior_face(int, const triangulation*, const recon_data*, double*);
void assemble_gradient_term_dirichlet_boundary_face(int, const triangulation*, const recon_data*, double*);
void evaluate_interior_face(int, const triangulation*, const recon_data*, const double*, double*);
void evaluate_boundary_face(int, const triangulation*, const recon_data*, const double*, double, double*);

/* Poisson equation solver. Solve the Poisson equation in the form -\nabla^2 u = f(x,y) */

enum poisson_solver_bcond { neumann, dirichlet };

typedef struct {
    triangulation* mesh;                // Mesh object associated with the solver
    double* phi;                        // Solution of the Poisson equation at cell centers of the mesh
    recon_data* clsq;                   // Constrained least-squares data structure on each face
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
void poisson_solver_construct_recon_data(poisson_solver*);
double poisson_solver_rhs(double,double);
double poisson_solver_exact_solution(double,double);
void poisson_solver_boundary_conditions(double, double, int, enum poisson_solver_bcond*,double*);
void poisson_solver_plot_vtk(poisson_solver*, const char*);
void poisson_solver_assemble_system(poisson_solver*);
void poisson_solver_run(poisson_solver*);

/* Unstready Incompressible Navier-Stokes Solver */

enum ins_solver_bcond { wall, inlet, outlet, symmetry };

typedef struct {
    double nu;                      // Kinematic viscosity of the fluid (constant)
    double CFL;                     // CFL constant (default 0.8)
    double dt;                      // Time step
    double time;                    // Time
    double tend;                    // Final time up to which the simulation needs to be done
    int timestep;                   // Time step number
    triangulation* mesh;            // Mesh object associated with the solver
    double *u, *v, *p;              // Velocity and pressure at cell centers
    double *u0, *v0;                // Velocities at previous time steps
    double *mdotf;                  // Mass flow rate through faces of the mesh
    double *uf, *vf, *pf;           // Velocity and pressure at faces of the mesh
    double *rhsu, *rhsv;            // RHS values of velocity for updating auxiliary momentum equation
    double **du, **dv;              // Gradient of velocity at cell centers
    double *divU;                   // Divergence of vellocity at cell centers of the mesh
    double *sp;                     // Source term for pressure poisson equation
    recon_data* clsq;               // Constrained least-squares reconstruction data for interior faces of the mesh (for all variables)
    recon_data* clsq_bnd_vel;       // Constrained least-squares reconstruction data for boundary faces of the mesh (for both u and v)
    recon_data* clsq_bnd_prs;       // Constrained least-squares reconstruction data for boundary faces of the mesh (for both pressure)
    enum ins_solver_bcond* bcond;   // Boundary condition associated with each boundary face
    double *ubnd, *vbnd, *pbnd;     // Velocity at pressure at boundary faces of the mesh
    double max_div_vel;             // Maximum divergence of velocity in the domain
    LIS_MATRIX A;                   // System matrix for solving pressure poisson equation (PPE)
    LIS_VECTOR rhs;                 // RHS vector for solving PPE
    LIS_VECTOR sol;                 // Solution vector of the PPE
    LIS_SOLVER solver;              // Iterative solver object
    int n_ppe_iter;                 // No. of iterations for PPE solver to converge
    double residual;                // Velocity residual (L2 norm)
} ins_solver;

ins_solver* ins_solver_allocate(const char*, double, double);
void ins_solver_initial_condition(double, double, double*, double*);
void ins_solver_boundary_conditions(double, double, int, enum ins_solver_bcond*, double*, double*, double*);
void ins_solver_compute_reconstruction_data(ins_solver*);
void ins_solver_calc_time_step(ins_solver*);
void ins_solver_calc_auxiliary_velocity(ins_solver*);
void ins_solver_calc_div_vel(ins_solver*);
void ins_solver_wall_viscous_flux(const ins_solver*, int, double, double, double, double, double*);
void ins_solver_assemble_ppe_matrix(ins_solver*);
void ins_solver_solve_ppe(ins_solver*);
void ins_solver_correct_velocity(ins_solver*);
void ins_solver_plot_vtk(ins_solver*);
void ins_solver_calc_residue(ins_solver*);
void ins_solver_calc_force_coeffs(ins_solver*, int, double, double, double*, double*);
void ins_solver_run(ins_solver*);
void ins_solver_free(ins_solver*);


#endif /* FV_H_ */
