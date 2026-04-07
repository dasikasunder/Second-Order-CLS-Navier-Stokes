/*
 * reconstruction.c
 *      Author: sunder
 */

#include "fv.h"

//---------------------------------------------------------------------
// Constructs constrained least squares data for an
// interior face of the mesh.
//---------------------------------------------------------------------

void clsq_construct_interior_face_data(int iface, const triangulation* mesh, recon_data* clsq) {

    int icell, idof, icell_local, c0, c1, signum;
    double x0, y0, dx1, dy1, dx2, dy2, dxi, dyi, dxi_sq, dyi_sq, dxi_dyi, dx, dy;

    // GSL matrices for solving the linear system

    gsl_matrix *A = gsl_matrix_calloc(ndof+2,ndof+2); // The extra 2 are the lagrange multipliers
    gsl_vector *b = gsl_vector_alloc(ndof+2);
    gsl_vector *sol = gsl_vector_alloc(ndof+2);
    gsl_permutation *p = gsl_permutation_alloc(ndof+2);

    x0 = mesh->face[iface].x[0]; y0 = mesh->face[iface].x[1];    // Center of the face

    c0 = mesh->face[iface].c[0]; c1 = mesh->face[iface].c[1];    // Cells adjacent to the face

    dx1 = mesh->cell[c0].x[0] - x0; dy1 = mesh->cell[c0].x[1] - y0;
    dx2 = mesh->cell[c1].x[0] - x0; dy2 = mesh->cell[c1].x[1] - y0;

    dxi = 0.; dyi = 0.; dxi_sq = 0.; dyi_sq = 0.; dxi_dyi = 0.;

    for (icell_local = 0; icell_local < mesh->face[iface].nvnb; ++icell_local) {

        icell = mesh->face[iface].vnb[icell_local];

        dx = mesh->cell[icell].x[0] - x0;
        dy = mesh->cell[icell].x[1] - y0;

        dxi += dx;
        dyi += dy;
        dxi_sq  += dx*dx;
        dyi_sq  += dy*dy;
        dxi_dyi += dx*dy;
    }

    gsl_matrix_set(A, 0, 0, 2.0*((double)(mesh->face[iface].nvnb)));
    gsl_matrix_set(A, 0, 1, 2.0*dxi);
    gsl_matrix_set(A, 0, 2, 2.0*dyi);
    gsl_matrix_set(A, 0, 3, 1.0);
    gsl_matrix_set(A, 0, 4, 1.0);

    gsl_matrix_set(A, 1, 0, 2.0*dxi);
    gsl_matrix_set(A, 1, 1, 2.0*dxi_sq);
    gsl_matrix_set(A, 1, 2, 2.0*dxi_dyi);
    gsl_matrix_set(A, 1, 3, dx1);
    gsl_matrix_set(A, 1, 4, dx2);

    gsl_matrix_set(A, 2, 0, 2.0*dyi);
    gsl_matrix_set(A, 2, 1, 2.0*dxi_dyi);
    gsl_matrix_set(A, 2, 2, 2.0*dyi_sq);
    gsl_matrix_set(A, 2, 3, dy1);
    gsl_matrix_set(A, 2, 4, dy2);

    gsl_matrix_set(A, 3, 0, 1.);
    gsl_matrix_set(A, 3, 1, dx1);
    gsl_matrix_set(A, 3, 2, dy1);

    gsl_matrix_set(A, 4, 0, 1.);
    gsl_matrix_set(A, 4, 1, dx2);
    gsl_matrix_set(A, 4, 2, dy2);

    gsl_linalg_LU_decomp(A, p, &signum);

    //---------------------------------------------------------------

    gsl_vector_set_zero(b);
    gsl_vector_set(b, 3, 1.0);

    gsl_linalg_LU_solve(A, p, b, sol);

    for (idof = 0; idof < ndof; ++idof)
        clsq->sol[0][idof] = gsl_vector_get(sol,idof);

    gsl_vector_set_zero(b);
    gsl_vector_set(b, 4, 1.0);

    gsl_linalg_LU_solve(A, p, b, sol);

    for (idof = 0; idof < ndof; ++idof)
        clsq->sol[1][idof] = gsl_vector_get(sol,idof);

    for (icell_local = 0; icell_local < mesh->face[iface].nvnb; ++icell_local) {

        icell = mesh->face[iface].vnb[icell_local];

        gsl_vector_set_zero(b);

        dx = mesh->cell[icell].x[0] - x0;
        dy = mesh->cell[icell].x[1] - y0;

        gsl_vector_set(b, 0, 2.0);
        gsl_vector_set(b, 1, 2.0*dx);
        gsl_vector_set(b, 2, 2.0*dy);

        gsl_linalg_LU_solve(A, p, b, sol);

        for (idof = 0; idof < ndof; ++idof)
            clsq->sol[icell_local+2][idof] = gsl_vector_get(sol,idof);
    }

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(sol);
    gsl_permutation_free(p);
}

//---------------------------------------------------------------------
// Constructs constrained least squares data for a face on dirichlet
// boundary
//---------------------------------------------------------------------

void clsq_construct_dirichlet_boundary_face_data(int iface, const triangulation* mesh, recon_data* clsq) {

    int icell, idof, icell_local, c0, signum;
    double x0, y0, dx1, dy1, dxi, dyi, dxi_sq, dyi_sq, dxi_dyi, dx, dy;

    gsl_matrix *A = gsl_matrix_calloc(ndof+2,ndof+2); // The extra 2 are the lagrange multipliers
    gsl_vector *b = gsl_vector_alloc(ndof+2);
    gsl_vector *sol = gsl_vector_alloc(ndof+2);
    gsl_permutation *p = gsl_permutation_alloc(ndof+2);

    c0 = mesh->face[iface].c[0];

    x0 = mesh->face[iface].x[0]; y0 = mesh->face[iface].x[1];

    dx1 = mesh->cell[c0].x[0] - x0; dy1 = mesh->cell[c0].x[1] - y0;

    dxi = 0.; dyi = 0.; dxi_sq = 0.; dyi_sq = 0.; dxi_dyi = 0.;

    for (icell_local = 0; icell_local < mesh->face[iface].nvnb; ++icell_local) {

        icell = mesh->face[iface].vnb[icell_local];

        dx = mesh->cell[icell].x[0] - x0;
        dy = mesh->cell[icell].x[1] - y0;

        dxi += dx;
        dyi += dy;
        dxi_sq  += dx*dx;
        dyi_sq  += dy*dy;
        dxi_dyi += dx*dy;
    }

    gsl_matrix_set(A, 0, 0, 2.0*((double)(mesh->face[iface].nvnb)));
    gsl_matrix_set(A, 0, 1, 2.0*dxi);
    gsl_matrix_set(A, 0, 2, 2.0*dyi);
    gsl_matrix_set(A, 0, 3, 1.0);
    gsl_matrix_set(A, 0, 4, 1.0);

    gsl_matrix_set(A, 1, 0, 2.0*dxi);
    gsl_matrix_set(A, 1, 1, 2.0*dxi_sq);
    gsl_matrix_set(A, 1, 2, 2.0*dxi_dyi);
    gsl_matrix_set(A, 1, 3, dx1);
    gsl_matrix_set(A, 1, 4, 0.);

    gsl_matrix_set(A, 2, 0, 2.0*dyi);
    gsl_matrix_set(A, 2, 1, 2.0*dxi_dyi);
    gsl_matrix_set(A, 2, 2, 2.0*dyi_sq);
    gsl_matrix_set(A, 2, 3, dy1);
    gsl_matrix_set(A, 2, 4, 0.);

    gsl_matrix_set(A, 3, 0, 1.);
    gsl_matrix_set(A, 3, 1, dx1);
    gsl_matrix_set(A, 3, 2, dy1);

    gsl_matrix_set(A, 4, 0, 1.);

    gsl_linalg_LU_decomp(A, p, &signum);

    //---------------------------------------------------------------

    gsl_vector_set_zero(b);
    gsl_vector_set(b, 3, 1.0);

    gsl_linalg_LU_solve(A, p, b, sol);

    for (idof = 0; idof < ndof; ++idof)
        clsq->sol[0][idof] = gsl_vector_get(sol,idof);

    gsl_vector_set_zero(b);
    gsl_vector_set(b, 4, 1.0);

    gsl_linalg_LU_solve(A, p, b, sol);

    for (idof = 0; idof < ndof; ++idof)
        clsq->sol[1][idof] = gsl_vector_get(sol,idof);

    for (icell_local = 0; icell_local < mesh->face[iface].nvnb; ++icell_local) {

        icell = mesh->face[iface].vnb[icell_local];

        gsl_vector_set_zero(b);

        dx = mesh->cell[icell].x[0] - x0;
        dy = mesh->cell[icell].x[1] - y0;

        gsl_vector_set(b, 0, 2.0);
        gsl_vector_set(b, 1, 2.0*dx);
        gsl_vector_set(b, 2, 2.0*dy);

        gsl_linalg_LU_solve(A, p, b, sol);

        for (idof = 0; idof < ndof; ++idof)
            clsq->sol[icell_local+2][idof] = gsl_vector_get(sol,idof);
    }

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(sol);
    gsl_permutation_free(p);
}

//---------------------------------------------------------------------
// Construct constrained least squares data for a face on neumann
// boundary
//---------------------------------------------------------------------

void clsq_construct_neumann_boundary_face_data(int iface, const triangulation* mesh, recon_data* clsq) {

    int icell, idof, icell_local, c0, signum;
    double x0, y0, dx1, dy1, dxi, dyi, dxi_sq, dyi_sq, dxi_dyi, dx, dy;
    double nx, ny;

    gsl_matrix *A = gsl_matrix_calloc(ndof+2,ndof+2); // The extra 2 are the lagrange multipliers
    gsl_vector *b = gsl_vector_alloc(ndof+2);
    gsl_vector *sol = gsl_vector_alloc(ndof+2);
    gsl_permutation *p = gsl_permutation_alloc(ndof+2);

    c0 = mesh->face[iface].c[0];

    x0 = mesh->face[iface].x[0]; y0 = mesh->face[iface].x[1];

    nx = mesh->face[iface].n[0]; ny = mesh->face[iface].n[1];

    dx1 = mesh->cell[c0].x[0] - x0; dy1 = mesh->cell[c0].x[1] - y0;

    dxi = 0.; dyi = 0.; dxi_sq = 0.; dyi_sq = 0.; dxi_dyi = 0.;

    for (icell_local = 0; icell_local < mesh->face[iface].nvnb; ++icell_local) {

        icell = mesh->face[iface].vnb[icell_local];

        dx = mesh->cell[icell].x[0] - x0;
        dy = mesh->cell[icell].x[1] - y0;

        dxi += dx;
        dyi += dy;
        dxi_sq  += dx*dx;
        dyi_sq  += dy*dy;
        dxi_dyi += dx*dy;
    }

    gsl_matrix_set(A, 0, 0, 2.0*((double)(mesh->face[iface].nvnb)));
    gsl_matrix_set(A, 0, 1, 2.0*dxi);
    gsl_matrix_set(A, 0, 2, 2.0*dyi);
    gsl_matrix_set(A, 0, 3, 1.0);
    gsl_matrix_set(A, 0, 4, 0.0);

    gsl_matrix_set(A, 1, 0, 2.0*dxi);
    gsl_matrix_set(A, 1, 1, 2.0*dxi_sq);
    gsl_matrix_set(A, 1, 2, 2.0*dxi_dyi);
    gsl_matrix_set(A, 1, 3, dx1);
    gsl_matrix_set(A, 1, 4, nx);

    gsl_matrix_set(A, 2, 0, 2.0*dyi);
    gsl_matrix_set(A, 2, 1, 2.0*dxi_dyi);
    gsl_matrix_set(A, 2, 2, 2.0*dyi_sq);
    gsl_matrix_set(A, 2, 3, dy1);
    gsl_matrix_set(A, 2, 4, ny);

    gsl_matrix_set(A, 3, 0, 1.);
    gsl_matrix_set(A, 3, 1, dx1);
    gsl_matrix_set(A, 3, 2, dy1);

    gsl_matrix_set(A, 4, 1, nx);
    gsl_matrix_set(A, 4, 2, ny);

    gsl_linalg_LU_decomp(A, p, &signum);

    //---------------------------------------------------------------

    gsl_vector_set_zero(b);
    gsl_vector_set(b, 3, 1.0);

    gsl_linalg_LU_solve(A, p, b, sol);

    for (idof = 0; idof < ndof; ++idof)
        clsq->sol[0][idof] = gsl_vector_get(sol,idof);

    gsl_vector_set_zero(b);
    gsl_vector_set(b, 4, 1.0);

    gsl_linalg_LU_solve(A, p, b, sol);

    for (idof = 0; idof < ndof; ++idof)
        clsq->sol[1][idof] = gsl_vector_get(sol,idof);

    for (icell_local = 0; icell_local < mesh->face[iface].nvnb; ++icell_local) {

        icell = mesh->face[iface].vnb[icell_local];

        gsl_vector_set_zero(b);

        dx = mesh->cell[icell].x[0] - x0;
        dy = mesh->cell[icell].x[1] - y0;

        gsl_vector_set(b, 0, 2.0);
        gsl_vector_set(b, 1, 2.0*dx);
        gsl_vector_set(b, 2, 2.0*dy);

        gsl_linalg_LU_solve(A, p, b, sol);

        for (idof = 0; idof < ndof; ++idof)
            clsq->sol[icell_local+2][idof] = gsl_vector_get(sol,idof);
    }

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(sol);
    gsl_permutation_free(p);
}

//---------------------------------------------------------------------
// Evaluate Constrained least squares data for a face on interior of
// the mesh
//---------------------------------------------------------------------

void evaluate_interior_face(int iface, const triangulation* mesh, const recon_data *clsq, const double *phi, double* value) {

    int icell, icell_local, idof, iowner, ineighbor;
    double phi_local;

    iowner    = mesh->face[iface].c[0];
    ineighbor = mesh->face[iface].c[1];

    phi_local = phi[iowner];

    for (idof = 0; idof < ndof; ++idof)
        value[idof] = phi_local*clsq->sol[0][idof];

    phi_local = phi[ineighbor];

    for (idof = 0; idof < ndof; ++idof)
        value[idof] += phi_local*clsq->sol[1][idof];

    for (icell_local = 0; icell_local < mesh->face[iface].nvnb; ++icell_local) {

        icell = mesh->face[iface].vnb[icell_local];

        phi_local = phi[icell];

        for (idof = 0; idof < ndof; ++idof)
            value[idof] += phi_local*clsq->sol[icell_local+2][idof];
    }
}

//---------------------------------------------------------------------
// Evaluate Constrained least squares data for a face on a boundary
// (Works for both dirichlet and Neumann boundary)
//---------------------------------------------------------------------


void evaluate_boundary_face(int iface, const triangulation* mesh, const recon_data *clsq, const double *phi, double phi_b, double* value) {

    int icell, icell_local, idof, iowner;
    double phi_local;

    iowner    = mesh->face[iface].c[0];

    phi_local = phi[iowner];

    for (idof = 0; idof < ndof; ++idof)
        value[idof] = phi_local*clsq->sol[0][idof];

    phi_local = phi_b;

    for (idof = 0; idof < ndof; ++idof)
        value[idof] += phi_local*clsq->sol[1][idof];

    for (icell_local = 0; icell_local < mesh->face[iface].nvnb; ++icell_local) {

        icell = mesh->face[iface].vnb[icell_local];

        phi_local = phi[icell];

        for (idof = 0; idof < ndof; ++idof)
            value[idof] += phi_local*clsq->sol[icell_local+2][idof];
    }
}

//---------------------------------------------------------------------
// Assemble gradient term for an intrior face
//---------------------------------------------------------------------

void assemble_gradient_term_interior_face(int iface, const triangulation* mesh,
                                          const recon_data* clsq, double* coeffs) {

    int icell_local;

    double nx   = mesh->face[iface].n[0];
    double ny   = mesh->face[iface].n[1];
    double area = mesh->face[iface].area;

    coeffs[0] = -area*(clsq->sol[0][1]*nx + clsq->sol[0][2]*ny); // Contribution to c0 cell
    coeffs[1] = -area*(clsq->sol[1][1]*nx + clsq->sol[1][2]*ny); // Contribution to c1 cell

    for (icell_local = 0; icell_local < mesh->face[iface].nvnb; ++icell_local) { // Contribution of vertex neighbouring cells
        coeffs[icell_local+2] = -area*(clsq->sol[icell_local+2][1]*nx + clsq->sol[icell_local+2][2]*ny);

    }
}

//---------------------------------------------------------------------
// Assemble gradient term for a face on dirichlet boundary
//---------------------------------------------------------------------

void assemble_gradient_term_dirichlet_boundary_face(int iface, const triangulation* mesh,
                                          const recon_data* clsq, double* coeffs) {

    int icell_local;

    double nx   = mesh->face[iface].n[0];
    double ny   = mesh->face[iface].n[1];
    double area = mesh->face[iface].area;

    coeffs[0] = -area*(clsq->sol[0][1]*nx + clsq->sol[0][2]*ny); // Contribution to c0 cell
    coeffs[1] =  area*(clsq->sol[1][1]*nx + clsq->sol[1][2]*ny); // Contribution of boundary condition - No minus sign here, so that
                                                                 // contribution is added to the RHS

    for (icell_local = 0; icell_local < mesh->face[iface].nvnb; ++icell_local) // Contribution of vertex neighbouring cells
        coeffs[icell_local+2] = -area*(clsq->sol[icell_local+2][1]*nx + clsq->sol[icell_local+2][2]*ny);
}
