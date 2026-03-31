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

int main() {

    const char mesh_file[] = "mesh.su2";

    poisson_solver* ps = poisson_solver_allocate(mesh_file);
    poisson_solver_run(ps);
    poisson_solver_plot_vtk(ps, "sol.vtk");
    poisson_solver_free(ps);

    return 0;
}

