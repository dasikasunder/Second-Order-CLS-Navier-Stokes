/*
 * main.c
 *      Author: sunder
 */

#include "fv.h"

double poisson_solver_exact_solution(double x, double y) {
    return sin(M_PI*x)*sin(M_PI*y);
}

double poisson_solver_rhs(double x, double y) {
    return 2.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}

void poisson_solver_boundary_conditions(double x, double y, int btag, enum poisson_solver_bcond* bcond, double* value) {
    if (btag == 230) {
        *bcond = dirichlet;
        *value = sin(M_PI*x)*sin(M_PI*y);
    }

    else if (btag == 231) {
        *bcond = neumann;
        *value = M_PI*sin(M_PI*x)*cos(M_PI*y);
    }

    else {
        printf("Error. Boundary condition for tag %d is not defined\n.", btag);
        exit(EXIT_FAILURE);
    }

}

int main() {

    const char mesh_file[] = "mesh.su2";

    poisson_solver* ps = poisson_solver_allocate(mesh_file);
    poisson_solver_run(ps);
    poisson_solver_plot_vtk(ps, "sol.vtk");
    poisson_solver_free(ps);

    return 0;
}

