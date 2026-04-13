/*
 * poisson_solver.c
 *      Author: sunder
 */

#include "fv.h"

/* Allocate memory for Poisson solver object */

poisson_solver* poisson_solver_allocate(const char* mesh_file) {

    int iface, btag, count;
    double xf, yf;

    poisson_solver* ps = (poisson_solver*)malloc(sizeof(poisson_solver));
    ps->mesh = triangulation_allocate(mesh_file);

    allocate_1d(ps->phi, ps->mesh->n_cell);
    allocate_1d(ps->clsq, ps->mesh->n_face);
    allocate_1d(ps->phi_bnd, ps->mesh->n_bnd_face);
    allocate_1d(ps->bcond, ps->mesh->n_bnd_face);
    allocate_1d(ps->bmap, ps->mesh->n_face);

    for (iface = 0; iface < ps->mesh->n_face; ++iface)
        allocate_2d(ps->clsq[iface].sol, ps->mesh->face[iface].nvnb+2, ndof);

    lis_vector_create(0,&ps->f); lis_vector_set_size(ps->f,0,ps->mesh->n_cell);
    lis_vector_create(0,&ps->u); lis_vector_set_size(ps->u,0,ps->mesh->n_cell);
    lis_matrix_create(0,&ps->A); lis_matrix_set_size(ps->A,0,ps->mesh->n_cell);

    lis_solver_create(&ps->solver);

    lis_solver_set_option("-i cg -p ilu", ps->solver);
    lis_solver_set_option("-tol 1.0e-10 -maxiter 10000", ps->solver);

    // Set boundary conditions and values

    count = 0;

    for (iface = ps->mesh->n_int_face; iface < ps->mesh->n_face; ++iface) {

        xf    = ps->mesh->face[iface].x[0];
        yf    = ps->mesh->face[iface].x[1];
        btag  = ps->mesh->face[iface].btag;

        poisson_solver_boundary_conditions(xf, yf, btag,  &ps->bcond[count], &ps->phi_bnd[count]);
        ps->bmap[iface] = count;
        count++;
    }

    return ps;
}

/* Release memory allocated to the Poisson solver */

void poisson_solver_free(poisson_solver* ps) {

    int iface;

    for (iface = 0; iface < ps->mesh->n_face; ++iface)
        free_2d(ps->clsq[iface].sol);

    free_1d(ps->phi);
    free_1d(ps->clsq);
    lis_vector_destroy(ps->u); lis_vector_destroy(ps->f);
    lis_matrix_destroy(ps->A);
    lis_solver_destroy(ps->solver);
    triangulation_free(ps->mesh);
}

/* Construct constrained least-squares data for each face in the mesh */

void poisson_solver_construct_recon_data(poisson_solver* ps) {

    int iface;

    // 1) Loop over interior faces and construct clsq data for interior faces

    for (iface = 0; iface < ps->mesh->n_int_face; ++iface)
        clsq_construct_interior_face_data(iface, ps->mesh, &ps->clsq[iface]);

    // 2) Loop over boundary faces and construct clsq data for boundary faces

    for (iface = ps->mesh->n_int_face; iface < ps->mesh->n_face; ++iface) {

        if (ps->bcond[ps->bmap[iface]] == dirichlet)
            clsq_construct_dirichlet_boundary_face_data(iface, ps->mesh, &ps->clsq[iface]);
        else
            clsq_construct_neumann_boundary_face_data(iface, ps->mesh, &ps->clsq[iface]);
    }
}

/* Assemble the LHS matrix for solving the poisson equation */

void poisson_solver_assemble_system(poisson_solver* ps) {

    int iface, iowner, ineighbor, icell, icell_local;
    double phi_b;
    double coeffs[25];

    // Loop over all interior faces

    for (iface = 0; iface < ps->mesh->n_int_face; ++iface) {

        assemble_gradient_term_interior_face(iface, ps->mesh, &ps->clsq[iface], coeffs);

        iowner = ps->mesh->face[iface].c[0];
        ineighbor = ps->mesh->face[iface].c[1];

        // Add contribution to owner cell

        lis_matrix_set_value(LIS_ADD_VALUE, iowner, iowner,    coeffs[0], ps->A);
        lis_matrix_set_value(LIS_ADD_VALUE, iowner, ineighbor, coeffs[1], ps->A);

        // Add contribution to ineighbor cell

        lis_matrix_set_value(LIS_ADD_VALUE, ineighbor, iowner,    -coeffs[0], ps->A);
        lis_matrix_set_value(LIS_ADD_VALUE, ineighbor, ineighbor, -coeffs[1], ps->A);

        // Add contribution of vertex neighbour cells to owner and neighbour cells

        for (icell_local = 0; icell_local < ps->mesh->face[iface].nvnb; ++icell_local) {

            icell = ps->mesh->face[iface].vnb[icell_local];

            lis_matrix_set_value(LIS_ADD_VALUE, iowner,    icell,  coeffs[icell_local+2], ps->A);
            lis_matrix_set_value(LIS_ADD_VALUE, ineighbor, icell, -coeffs[icell_local+2], ps->A);
        }
    }

    // Loop over all boundary faces

    for (iface = ps->mesh->n_int_face; iface < ps->mesh->n_face; ++iface) {

        iowner = ps->mesh->face[iface].c[0];

        if (ps->bcond[ps->bmap[iface]] == dirichlet) {

            assemble_gradient_term_dirichlet_boundary_face(iface, ps->mesh, &ps->clsq[iface], coeffs);

            lis_matrix_set_value(LIS_ADD_VALUE, iowner, iowner, coeffs[0], ps->A);

            for (icell_local = 0; icell_local < ps->mesh->face[iface].nvnb; ++icell_local) {
                icell = ps->mesh->face[iface].vnb[icell_local];
                lis_matrix_set_value(LIS_ADD_VALUE, iowner,    icell,  coeffs[icell_local+2], ps->A);
            }

            phi_b = ps->phi_bnd[ps->bmap[iface]];
            lis_vector_set_value(LIS_ADD_VALUE, iowner, coeffs[1]*phi_b, ps->f);
        }

        else {
            phi_b = ps->phi_bnd[ps->bmap[iface]];
            lis_vector_set_value(LIS_ADD_VALUE, iowner, ps->mesh->face[iface].area*phi_b, ps->f);
        }
    }

    lis_matrix_set_type(ps->A,LIS_MATRIX_CSR);
    lis_matrix_assemble(ps->A);
}

/* Plot solution in VTK ASCII format */

void poisson_solver_plot_vtk(poisson_solver* ps, const char* out_file_name) {

    int ivrtx;
    double value;

    triangulation_plot_vtk(ps->mesh, out_file_name); // Plot mesh connectivity

    FILE *file_ptr;

    file_ptr = fopen(out_file_name, "a");

    if (file_ptr == NULL)
        printf("Error opening the file.\n");

    fprintf(file_ptr, "CELL_DATA %d\n", ps->mesh->n_cell);

    fprintf(file_ptr, "SCALARS Phi double 1\n");
    fprintf(file_ptr, "LOOKUP_TABLE default\n");

    for (ivrtx = 0; ivrtx < ps->mesh->n_cell; ++ivrtx) {
        lis_vector_get_value(ps->u, ivrtx, &value);
        fprintf(file_ptr, "%.5e\n", value);
    }

    fclose(file_ptr);
}

/* Put everything together and run the problem */

void poisson_solver_run(poisson_solver* ps) {

    int icell, n_iter;
    double xc, yc, vol;

    poisson_solver_construct_recon_data(ps); // Calculate CLSQ data for all the faces

    poisson_solver_assemble_system(ps);     // Assemble the Poisson system

    // Construct RHS

    for (icell = 0; icell < ps->mesh->n_cell; ++icell) {

        vol = ps->mesh->cell[icell].vol;
        xc = ps->mesh->cell[icell].x[0];
        yc = ps->mesh->cell[icell].x[1];

        lis_vector_set_value(LIS_ADD_VALUE, icell, vol*poisson_solver_rhs(xc,yc), ps->f);

    } // Cell loop

    // Solve system

    lis_solve(ps->A,ps->f,ps->u,ps->solver);
    lis_solver_get_iter(ps->solver, &n_iter);

    printf("Solver has converged in %d iterations\n", n_iter);

    // Calculate error

    double max_error = 0.0;
    double l2_error = 0.0;
    double value, error;

    for (icell = 0; icell < ps->mesh->n_cell; ++icell) {

        xc = ps->mesh->cell[icell].x[0];
        yc = ps->mesh->cell[icell].x[1];

        double exact = poisson_solver_exact_solution(xc,yc);
        lis_vector_get_value(ps->u, icell, &value);
        error = fabs(exact - value);

        l2_error += error*error;

        if (error > max_error)
            max_error = error;
    }

    l2_error = sqrt(l2_error/(double)ps->mesh->n_cell);

    printf("L2 Error  = %.3e\n", l2_error);
    printf("Max Error = %.3e\n", max_error);
}

double poisson_solver_exact_solution(double x, double y) {

    double R = 1.0;
    double r = hypot(x,y);
    double theta = atan2(y,x);
    double U = 1.0;

    return U*(r+R*R/r)*cos(theta);
}

double poisson_solver_rhs(double x, double y) {
    return 0.0*x+0.0*y;
}

void poisson_solver_boundary_conditions(double x, double y, int btag, enum poisson_solver_bcond* bcond, double* value) {

    if (btag == 4) {
        *bcond = dirichlet;
        *value = poisson_solver_exact_solution(x,y);
    }

    else if (btag == 19) {
        *bcond = neumann;
        *value = 0.0;
    }

    else {
        printf("Error. Boundary condition for tag %d is not defined\n.", btag);
        exit(EXIT_FAILURE);
    }

}
