/*
 * main.c
 *      Author: sunder
 */

#include "fv.h"

void ins_solver_initial_condition(double x, double y, double *u, double *v) {
    (void) x; (void) y; // Just to avoid unused variable warning
    *u = 0.; *v = 0.;
}

/*
void ins_solver_boundary_conditions(double x, double y, int btag, enum ins_solver_bcond* bc,
                                    double* u, double* v, double* p) {

    (void) x; (void) y;

    switch (btag) {
        case 230:
            *bc = wall;
            *u = 0.0; *v = 0.0; *p = 0.0;
            break;
        case 231:
            *bc = wall;
            *u = 1.0; *v = 0.0; *p = 0.0;
            break;
        default:
            printf("Error. Boundary condition %d not present\n", btag);
            exit(EXIT_FAILURE);
    }
}
*/

void ins_solver_boundary_conditions(double x, double y, int btag, enum ins_solver_bcond* bc,
                                    double* u, double* v, double* p) {

    (void) x; (void) y;

    switch (btag) {
        case 23:
            *bc = wall;
            *u = 0.0; *v = 0.0; *p = 0.0;
            break;
        case 9:
            *bc = inlet;
            *u = 24.*y*(0.5-y); *v = 0.0; *p = 0.0;
            *u = 4.0*0.3*y*(0.41-y)/(0.41*0.41);
            break;
        case 15:
            *bc = outlet;
            *u = 0.0; *v = 0.0; *p = 0.0;
            break;
        default:
            printf("Error. Boundary condition %d not present\n", btag);
            exit(EXIT_FAILURE);
    }
}


int main() {


    const char mesh_file[] = "mesh.su2"; // Name of the mesh file
    double nu = 1.0e-2;                 // Kinematic viscosity of the fluid
    double tend = 20.0;                  // Final time

    ins_solver* ins = ins_solver_allocate(mesh_file,nu,tend);
    ins_solver_run(ins);
    ins_solver_free(ins);

    /*
    poisson_solver* ps = poisson_solver_allocate(mesh_file);
    poisson_solver_run(ps);
    poisson_solver_plot_vtk(ps, "sol.vtk");
    poisson_solver_free(ps);
    */


    return 0;
}

