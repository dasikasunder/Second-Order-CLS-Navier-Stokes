/*
 * ins_solver.c
 *      Author: sunder
 */

#include "fv.h"

//---------------------------------------------------------------------
// Allocate memory for incompressible Navier-Stokes equation solver.
// Initialize the solution with the given initial conditions
//---------------------------------------------------------------------

ins_solver* ins_solver_allocate(const char* mesh_file, double nu, double tend) {

    int icell, iface, btag, bcount;
    double xf, yf;

    ins_solver* ins = (ins_solver*)malloc(sizeof(ins_solver)); // Create the ins_solver object
    ins->nu = nu;                                              // Set kinematic viscosity of the fluid (constant)
    ins->CFL = 0.6;                                            // CFL constant
    ins->dt = 0.0;                                             // Time step size (is calculated based on CFL condition and solution)
    ins->time = 0.0;                                           // Simulation time
    ins->tend = tend;                                          // Final time of simulation
    ins->timestep = 0;                                         // Time step number
    ins->mesh = triangulation_allocate(mesh_file);             // Initialize the mesh

    // Allocate memory

    allocate_1d(ins->u, ins->mesh->n_cell);                // x-velocity at cell centers
    allocate_1d(ins->v, ins->mesh->n_cell);                // y-velocity at cell centers
    allocate_1d(ins->p, ins->mesh->n_cell);                // Pressure at cell centers
    allocate_1d(ins->u0, ins->mesh->n_cell);               // x-velocity at cell centers (previous time step)
    allocate_1d(ins->v0, ins->mesh->n_cell);               // y-velocity at cell centers (previous time step)
    allocate_1d(ins->mdotf, ins->mesh->n_face);            // Mass flow rate through faces of the mesh
    allocate_1d(ins->uf, ins->mesh->n_face);               // x-velocity at face centers
    allocate_1d(ins->vf, ins->mesh->n_face);               // y-velocity at face centers
    allocate_1d(ins->pf, ins->mesh->n_face);               // Pressure at face centers
    allocate_1d(ins->rhsu, ins->mesh->n_cell);             // RHS value of x-velocity for updating auxiliary momentum equation
    allocate_1d(ins->rhsv, ins->mesh->n_cell);             // RHS value of y-velocity for updating auxiliary momentum equation
    allocate_2d(ins->du, ins->mesh->n_cell, 2);            // Gradient of x-velocity at cell centers
    allocate_2d(ins->dv, ins->mesh->n_cell, 2);            // Gradient of y-velocity at cell centers
    allocate_1d(ins->divU, ins->mesh->n_cell);             // Divergence of velocity at cell centers of the mesh
    allocate_1d(ins->sp, ins->mesh->n_cell);               // Source term for pressure poisson equation
    allocate_1d(ins->clsq,  ins->mesh->n_int_face);        // Reconstruction data on interior faces
    allocate_1d(ins->clsq_bnd_vel, ins->mesh->n_bnd_face); // Velocity reconstruction data on boundary faces
    allocate_1d(ins->clsq_bnd_prs, ins->mesh->n_bnd_face); // Pressure reconstruction data on boundary faces
    allocate_1d(ins->bcond, ins->mesh->n_bnd_face);        // Boundary condition associated with boundary faces
    allocate_1d(ins->ubnd, ins->mesh->n_bnd_face);         // x-velocity at boundary faces
    allocate_1d(ins->vbnd, ins->mesh->n_bnd_face);         // y-velocity at boundary faces
    allocate_1d(ins->pbnd, ins->mesh->n_bnd_face);         // Pressure at boundary faces

    // Initialize linear solver related data

    lis_matrix_create(0,&ins->A);   lis_matrix_set_size(ins->A,  0,ins->mesh->n_cell);
    lis_vector_create(0,&ins->rhs); lis_vector_set_size(ins->rhs,0,ins->mesh->n_cell);
    lis_vector_create(0,&ins->sol); lis_vector_set_size(ins->sol,0,ins->mesh->n_cell);

    lis_solver_create(&ins->solver);

    lis_solver_set_option("-i cg -p ilu", ins->solver);
    lis_solver_set_option("-tol 1.0e-8 -maxiter 10000", ins->solver);

    lis_vector_set_all(0.0,ins->sol);

    for (iface = 0; iface < ins->mesh->n_int_face; ++iface)
        allocate_2d(ins->clsq[iface].sol, ins->mesh->face[iface].nvnb+2, ndof);

    // Initialize solution at cell centers

    for (icell = 0; icell < ins->mesh->n_cell; ++icell) {
        ins_solver_initial_condition(ins->mesh->cell[icell].x[0],     // x-coordinate of cell center
                                     ins->mesh->cell[icell].x[1],     // y-coordinate of cell center
                                     &ins->u[icell], &ins->v[icell]); // x and y velocities
        ins->p[icell] = 0.0;                                          // pressure
        ins->sp[icell] = 0.0;

    }

    // Enumerate boundary conditions

    bcount = 0;

    for (iface = ins->mesh->n_int_face; iface < ins->mesh->n_face; ++iface) {

        allocate_2d(ins->clsq_bnd_vel[bcount].sol, ins->mesh->face[iface].nvnb+2, ndof);
        allocate_2d(ins->clsq_bnd_prs[bcount].sol, ins->mesh->face[iface].nvnb+2, ndof);

        xf    = ins->mesh->face[iface].x[0];
        yf    = ins->mesh->face[iface].x[1];
        btag  = ins->mesh->face[iface].btag;

        ins_solver_boundary_conditions(xf, yf, btag,  &ins->bcond[bcount],
                                       &ins->ubnd[bcount],
                                       &ins->vbnd[bcount],
                                       &ins->pbnd[bcount]);
        bcount++;
    }

    // Initialize solution face centers

    for (iface = 0; iface < ins->mesh->n_int_face; ++iface) {
        ins->uf[iface] = 0.0; ins->vf[iface] = 0.0;
        ins->pf[iface] = 0.0; ins->mdotf[iface] = 0.0;
    }

    bcount = 0;

    for (iface = ins->mesh->n_int_face; iface < ins->mesh->n_face; ++iface) {

        ins->uf[iface] = ins->ubnd[bcount]; ins->vf[iface] = ins->vbnd[bcount];
        ins->pf[iface] = ins->pbnd[bcount]; ins->mdotf[iface] = 0.0;

        bcount++;
    }

    return ins;
}

//---------------------------------------------------------------------
// Release all the memory allocated to the incompressible
// Navier-Stokes solver
//---------------------------------------------------------------------

void ins_solver_free(ins_solver* ins) {

    int iface;

    free_1d(ins->u); free_1d(ins->v); free_1d(ins->p);
    free_1d(ins->u0); free_1d(ins->v0);
    free_1d(ins->mdotf); free_1d(ins->rhsu); free_1d(ins->rhsv);
    free_2d(ins->du); free_2d(ins->dv);
    free_1d(ins->divU); free_1d(ins->sp);

    for (iface = 0; iface < ins->mesh->n_int_face; ++iface)
        free_2d(ins->clsq[iface].sol);

    for (iface = ins->mesh->n_int_face; iface < ins->mesh->n_face; ++iface) {
        free_2d(ins->clsq_bnd_vel[iface-ins->mesh->n_int_face].sol);
        free_2d(ins->clsq_bnd_prs[iface-ins->mesh->n_int_face].sol);
    }

    free_1d(ins->clsq); free_1d(ins->clsq_bnd_vel); free_1d(ins->clsq_bnd_prs);
    free_1d(ins->bcond); free_1d(ins->ubnd); free_1d(ins->vbnd); free_1d(ins->pbnd);
    free_1d(ins->uf); free_1d(ins->vf); free_1d(ins->pf);
    triangulation_free(ins->mesh);
    lis_vector_destroy(ins->sol); lis_vector_destroy(ins->rhs);
    lis_matrix_destroy(ins->A);
    lis_solver_destroy(ins->solver);
    free(ins);
}

//---------------------------------------------------------------------
// Precompute reconstruction data for both interior faces and boundary
// faces
//---------------------------------------------------------------------

void ins_solver_compute_reconstruction_data(ins_solver* ins) {

    int iface, counter;

    // 1) Loop over interior faces and construct clsq data for interior faces

    for (iface = 0; iface < ins->mesh->n_int_face; ++iface)
        clsq_construct_interior_face_data(iface, ins->mesh, &ins->clsq[iface]);

    // 2) Loop over boundary faces and construct clsq data according to the boundary condition

    counter = 0;

    for (iface = ins->mesh->n_int_face; iface < ins->mesh->n_face; ++iface) {

        if (ins->bcond[counter] == wall || ins->bcond[counter] == inlet) { // Velocity -> Dirichlet, Pressure -> Neumann
            clsq_construct_dirichlet_boundary_face_data(iface, ins->mesh, &ins->clsq_bnd_vel[counter]);
            clsq_construct_neumann_boundary_face_data(iface, ins->mesh, &ins->clsq_bnd_prs[counter]);
        }

        else if (ins->bcond[counter] == outlet) { // Velocity -> Neumann, Pressure -> Dirchlet
            clsq_construct_neumann_boundary_face_data(iface, ins->mesh, &ins->clsq_bnd_vel[counter]);
            clsq_construct_dirichlet_boundary_face_data(iface, ins->mesh, &ins->clsq_bnd_prs[counter]);
        }

        else {
            printf("Error. Specified boundary condition is not implemented \n");
            exit(EXIT_FAILURE);
        }
        counter++;
    }
}

//---------------------------------------------------------------------
// Calculate time step according to CFL condition
//---------------------------------------------------------------------

void ins_solver_calc_time_step(ins_solver* ins) {

    int icell;
    double h, vel, dt_local;

    ins->dt = 1.0e5; // Set arbitrarily large number

    for (icell = 0; icell < ins->mesh->n_cell; ++icell) {

        h = sqrt(ins->mesh->cell[icell].vol);            // Effective cell space based on sqrt(volume)
        vel = hypot(ins->u[icell], ins->v[icell]);       // Magnitude of velocity on the face
        dt_local = 0.5*ins->CFL*h/(vel + 2.0*ins->nu/h); // Local time step on the face (0.5 is for 2D)
        if (dt_local < ins->dt)
            ins->dt = dt_local;
    }
}

//---------------------------------------------------------------------
// Assemble system for solving pressure poisson equation
//---------------------------------------------------------------------

void ins_solver_assemble_ppe_matrix(ins_solver* ins) {

    int iface, iowner, ineighbor, icell, icell_local, bcount;
    double pb;
    double coeffs[25];

    // Loop over all interior faces

    for (iface = 0; iface < ins->mesh->n_int_face; ++iface) {

        assemble_gradient_term_interior_face(iface, ins->mesh, &ins->clsq[iface], coeffs);

        iowner = ins->mesh->face[iface].c[0];
        ineighbor = ins->mesh->face[iface].c[1];

        // Add contribution to owner cell

        lis_matrix_set_value(LIS_ADD_VALUE, iowner, iowner,    coeffs[0], ins->A);
        lis_matrix_set_value(LIS_ADD_VALUE, iowner, ineighbor, coeffs[1], ins->A);

        // Add contribution to ineighbor cell

        lis_matrix_set_value(LIS_ADD_VALUE, ineighbor, iowner,    -coeffs[0], ins->A);
        lis_matrix_set_value(LIS_ADD_VALUE, ineighbor, ineighbor, -coeffs[1], ins->A);

        // Add contribution of vertex neighbour cells to owner and neighbour cells

        for (icell_local = 0; icell_local < ins->mesh->face[iface].nvnb; ++icell_local) {

            icell = ins->mesh->face[iface].vnb[icell_local];

            lis_matrix_set_value(LIS_ADD_VALUE, iowner,    icell,  coeffs[icell_local+2], ins->A);
            lis_matrix_set_value(LIS_ADD_VALUE, ineighbor, icell, -coeffs[icell_local+2], ins->A);
        }
    }

    // Loop over all boundary faces

    bcount = 0;

    for (iface = ins->mesh->n_int_face; iface < ins->mesh->n_face; ++iface) {

        iowner = ins->mesh->face[iface].c[0];

        if (ins->bcond[bcount] == outlet) {

            assemble_gradient_term_dirichlet_boundary_face(iface, ins->mesh, &ins->clsq_bnd_prs[bcount], coeffs);

            lis_matrix_set_value(LIS_ADD_VALUE, iowner, iowner, coeffs[0], ins->A);

            for (icell_local = 0; icell_local < ins->mesh->face[iface].nvnb; ++icell_local) {
                icell = ins->mesh->face[iface].vnb[icell_local];
                lis_matrix_set_value(LIS_ADD_VALUE, iowner,    icell,  coeffs[icell_local+2], ins->A);
            }

            pb = ins->pbnd[bcount];
            ins->sp[iowner] += coeffs[1]*pb;
        }

        else {
            pb = ins->pbnd[bcount];
            ins->sp[iowner] += ins->mesh->face[iface].area*pb;
        }

        bcount++;
    }

    lis_matrix_set_type(ins->A,LIS_MATRIX_CSR);
    lis_matrix_assemble(ins->A);
}

//---------------------------------------------------------------------
// Calculate auxiliary velocity (velocity field without considering)
// pressure in the momentum equation)
//---------------------------------------------------------------------

void ins_solver_calc_auxiliary_velocity(ins_solver* ins) {

    int icell, iface, c0, c1, bcount;
    double dudx, dudy, dvdx, dvdy, flux, ul, ur, vl, vr;
    double dx, dy, xf, yf, nx, ny, area, vol, trm, value[ndof];
    double fv[2];

    // 0) Reset a number of values to zero

    for (icell = 0; icell < ins->mesh->n_cell; ++icell) {
        ins->rhsu[icell] = 0.0; ins->rhsv[icell] = 0.0;
        ins->du[icell][0] = 0.0; ins->du[icell][1] = 0.0;
        ins->dv[icell][0] = 0.0; ins->dv[icell][1] = 0.0;
        ins->u0[icell] = ins->u[icell]; ins->v0[icell] = ins->v[icell];
        ins->divU[icell] = 0.0;
    }

    // 1) Calculate viscous fluxes

    // 1a) -> Interior faces

    for (iface = 0; iface < ins->mesh->n_int_face; ++iface) {

        area = ins->mesh->face[iface].area;
        nx   = ins->mesh->face[iface].n[0]; ny   = ins->mesh->face[iface].n[1];
        c0   = ins->mesh->face[iface].c[0]; c1   = ins->mesh->face[iface].c[1];

        evaluate_interior_face(iface, ins->mesh, &ins->clsq[iface], ins->u, value);
        ins->uf[iface] = value[0]; dudx = value[1]; dudy = value[2];

        fv[0] = -ins->nu*(dudx*nx + dudy*ny)*area;

        ins->rhsu[c0] += fv[0]; ins->rhsu[c1] -= fv[0];

        evaluate_interior_face(iface, ins->mesh, &ins->clsq[iface], ins->v, value);
        ins->vf[iface] = value[0]; dvdx = value[1]; dvdy = value[2];

        fv[1] = -ins->nu*(dvdx*nx + dvdy*ny)*area;

        ins->rhsv[c0] += fv[1]; ins->rhsv[c1] -= fv[1];
    }

    // 1b) -> Boundary faces

    bcount = 0;

    for (iface = ins->mesh->n_int_face; iface < ins->mesh->n_face; ++iface) {

        area = ins->mesh->face[iface].area;
        c0   = ins->mesh->face[iface].c[0];
        nx   = ins->mesh->face[iface].n[0]; ny   = ins->mesh->face[iface].n[1];


        evaluate_boundary_face(iface, ins->mesh, &ins->clsq_bnd_vel[bcount], ins->u, ins->ubnd[bcount], value);
        ins->uf[iface] = value[0]; dudx = value[1]; dudy = value[2];

        evaluate_boundary_face(iface, ins->mesh, &ins->clsq_bnd_vel[bcount], ins->v, ins->vbnd[bcount], value);
        ins->vf[iface] = value[0]; dvdx = value[1]; dvdy = value[2];


        if (ins->bcond[bcount] == wall)
            ins_solver_wall_viscous_flux(ins, iface, dudx, dudy, dvdx, dvdy, fv);
        else {
            fv[0] = -ins->nu*(dudx*nx + dudy*ny)*area;
            fv[1] = -ins->nu*(dvdx*nx + dvdy*ny)*area;
        }

        ins->rhsu[c0] += fv[0];
        ins->rhsv[c0] += fv[1];

        bcount++;
    }

    // 2) Calculate convective fluxes

    // 2a) First calculate velocity gradients at cell centers using Green-Gauss theorem

    for (iface = 0; iface < ins->mesh->n_face; ++iface) {

        area = ins->mesh->face[iface].area;
        nx   = ins->mesh->face[iface].n[0]; ny = ins->mesh->face[iface].n[1];
        c0   = ins->mesh->face[iface].c[0]; c1 = ins->mesh->face[iface].c[1];

        trm = ins->uf[iface]*area;

        ins->du[c0][0] += trm*nx/ins->mesh->cell[c0].vol;
        ins->du[c0][1] += trm*ny/ins->mesh->cell[c0].vol;

        if (! ins->mesh->face[iface].at_boundary) {
            ins->du[c1][0] -= trm*nx/ins->mesh->cell[c1].vol;
            ins->du[c1][1] -= trm*ny/ins->mesh->cell[c1].vol;
        }

        trm = ins->vf[iface]*area;

        ins->dv[c0][0] += trm*nx/ins->mesh->cell[c0].vol;
        ins->dv[c0][1] += trm*ny/ins->mesh->cell[c0].vol;

        if (! ins->mesh->face[iface].at_boundary) {
            ins->dv[c1][0] -= trm*nx/ins->mesh->cell[c1].vol;
            ins->dv[c1][1] -= trm*ny/ins->mesh->cell[c1].vol;
        }
    }

    // 2b) Calculate convective flux through interior faces of the mesh

    for (iface = 0; iface < ins->mesh->n_int_face; ++iface) {

        xf = ins->mesh->face[iface].x[0]; yf = ins->mesh->face[iface].x[1];
        c0   = ins->mesh->face[iface].c[0]; c1 = ins->mesh->face[iface].c[1];

        dx = xf - ins->mesh->cell[c0].x[0]; dy = yf - ins->mesh->cell[c0].x[1];

        ul = ins->u[c0] + dx*ins->du[c0][0] + dy*ins->du[c0][1];
        vl = ins->v[c0] + dx*ins->dv[c0][0] + dy*ins->dv[c0][1];

        dx = xf - ins->mesh->cell[c1].x[0]; dy = yf - ins->mesh->cell[c1].x[1];

        ur = ins->u[c1] + dx*ins->du[c1][0] + dy*ins->du[c1][1];
        vr = ins->v[c1] + dx*ins->dv[c1][0] + dy*ins->dv[c1][1];

        if (ins->mdotf[iface] > 0.0) {
            flux = ins->mdotf[iface]*ul;
            ins->rhsu[c0] += flux; ins->rhsu[c1] -= flux;
            flux = ins->mdotf[iface]*vl;
            ins->rhsv[c0] += flux; ins->rhsv[c1] -= flux;
        }

        else {
            flux = ins->mdotf[iface]*ur;
            ins->rhsu[c0] += flux; ins->rhsu[c1] -= flux;
            flux = ins->mdotf[iface]*vr;
            ins->rhsv[c0] += flux; ins->rhsv[c1] -= flux;
        }
    }

    // 2c) Calculate convective flux through boundary faces

    bcount = 0;

    for (iface = ins->mesh->n_int_face; iface < ins->mesh->n_face; ++iface) {

        xf = ins->mesh->face[iface].x[0]; yf = ins->mesh->face[iface].x[1];
        c0   = ins->mesh->face[iface].c[0];

        if (ins->bcond[bcount] == wall || ins->bcond[bcount] == inlet) {
            ul = ins->ubnd[bcount]; vl = ins->vbnd[bcount];
        }

        else {

            dx = xf - ins->mesh->cell[c0].x[0]; dy = yf - ins->mesh->cell[c0].x[1];

            ul = ins->u[c0] + dx*ins->du[c0][0] + dy*ins->du[c0][1];
            vl = ins->v[c0] + dx*ins->dv[c0][0] + dy*ins->dv[c0][1];
        }

        flux = ins->mdotf[iface]*ul; ins->rhsu[c0] += flux;
        flux = ins->mdotf[iface]*vl; ins->rhsv[c0] += flux;

        bcount++;
    }

    // 3) Update cell center velocity

    for (icell = 0; icell < ins->mesh->n_cell; ++icell) {
        vol = ins->mesh->cell[icell].vol;
        ins->u[icell] -= ins->dt*ins->rhsu[icell]/vol;
        ins->v[icell] -= ins->dt*ins->rhsv[icell]/vol;
    }
}

//---------------------------------------------------------------------
// Viscous flux at a wall
//---------------------------------------------------------------------

void ins_solver_wall_viscous_flux(const ins_solver* ins, int iface, double dudx, double dudy, double dvdx, double dvdy, double* flux) {

    double nu = ins->nu;
    double area = ins->mesh->face[iface].area;
    double nx = ins->mesh->face[iface].n[0]; double ny = ins->mesh->face[iface].n[1];

    double tauxx = 2.0*nu*dudx;
    double tauxy = nu*(dudy + dvdx);
    double tauyy = 2.0*nu*dvdy;

    flux[0] = -area*(tauxx*nx + tauxy*ny);
    flux[1] = -area*(tauxy*nx + tauyy*ny);
}

//---------------------------------------------------------------------
// Calculate divergence of velocity in each cell using Gauss-Divergence
// theorem
//---------------------------------------------------------------------

void ins_solver_calc_div_vel(ins_solver* ins) {

    int iface, c0, c1, bcount;
    double mdotf, area, u, v,  nx, ny;
    double value[3];

    for (iface = 0; iface < ins->mesh->n_int_face; ++iface) {

        area = ins->mesh->face[iface].area;
        nx   = ins->mesh->face[iface].n[0]; ny   = ins->mesh->face[iface].n[1];
        c0   = ins->mesh->face[iface].c[0]; c1   = ins->mesh->face[iface].c[1];

        evaluate_interior_face(iface, ins->mesh, &ins->clsq[iface], ins->u, value); u = value[0];
        evaluate_interior_face(iface, ins->mesh, &ins->clsq[iface], ins->v, value); v = value[0];

        ins->uf[iface] = u; ins->vf[iface] = v;

        mdotf = (u*nx + v*ny)*area;

        ins->divU[c0] += mdotf/ins->mesh->cell[c0].vol;
        ins->divU[c1] -= mdotf/ins->mesh->cell[c1].vol;
    }

    bcount = 0;

    for (iface = ins->mesh->n_int_face; iface < ins->mesh->n_face; ++iface) {

        area = ins->mesh->face[iface].area;
        nx   = ins->mesh->face[iface].n[0]; ny   = ins->mesh->face[iface].n[1];
        c0   = ins->mesh->face[iface].c[0];

        if (ins->bcond[bcount] == wall)
            mdotf = 0.0;

        else if (ins->bcond[bcount] == inlet)
            mdotf = (ins->ubnd[bcount]*nx + ins->vbnd[bcount]*ny)*area;

        else {

            evaluate_boundary_face(iface, ins->mesh, &ins->clsq_bnd_vel[bcount], ins->u, ins->ubnd[bcount], value); u = value[0];
            evaluate_boundary_face(iface, ins->mesh, &ins->clsq_bnd_vel[bcount], ins->v, ins->vbnd[bcount], value); v = value[0];

            ins->uf[iface] = u; ins->vf[iface] = v;

            mdotf = (u*nx + v*ny)*area;
        }

        ins->divU[c0] += mdotf/ins->mesh->cell[c0].vol;

        bcount++;
    }
}

//---------------------------------------------------------------------
// Solve pressure poisson equation
//---------------------------------------------------------------------

void ins_solver_solve_ppe(ins_solver* ins) {

    int icell;
    double value, vol;

    // Set RHS for the matrix

    for (icell = 0; icell < ins->mesh->n_cell; ++icell) {

        vol = ins->mesh->cell[icell].vol;
        value = ins->sp[icell] - vol*ins->divU[icell]/ins->dt;

        lis_vector_set_value(LIS_INS_VALUE, icell, value, ins->rhs);
    }

    lis_solve(ins->A,ins->rhs,ins->sol,ins->solver);
    lis_solver_get_iter(ins->solver, &ins->n_ppe_iter);

    // Store solution in p array

    for (icell = 0; icell < ins->mesh->n_cell; ++icell) {
        lis_vector_get_value(ins->sol, icell, &value);
        ins->p[icell] = value;
    }
}

//---------------------------------------------------------------------
// Correct velocity at cell centers and face centers
//---------------------------------------------------------------------

void ins_solver_correct_velocity(ins_solver* ins) {

    int iface, icell, bcount;
    double area, nx, ny, un, dpdn, value[3];
    double div_vel, vol, snsign, dp[2];

    // 1a) Correct velocity at interior faces of the mesh

    for (iface = 0; iface < ins->mesh->n_int_face; ++iface) {

        area = ins->mesh->face[iface].area;
        nx   = ins->mesh->face[iface].n[0]; ny   = ins->mesh->face[iface].n[1];

        evaluate_interior_face(iface, ins->mesh, &ins->clsq[iface], ins->p, value);

        ins->pf[iface] = value[0];

        dpdn = value[1]*nx + value[2]*ny;

        un = ins->uf[iface]*nx + ins->vf[iface]*ny;

        ins->mdotf[iface] = (un - ins->dt*dpdn)*area;
    }

    // 1b) Correct velocity at boundary faces of the mesh

    bcount = 0;

    for (iface = ins->mesh->n_int_face; iface < ins->mesh->n_face; ++iface) {

        area = ins->mesh->face[iface].area;
        nx   = ins->mesh->face[iface].n[0]; ny   = ins->mesh->face[iface].n[1];

        evaluate_boundary_face(iface, ins->mesh, &ins->clsq_bnd_prs[bcount], ins->p, ins->pbnd[bcount], value);

        ins->pf[iface] = value[0];

        if (ins->bcond[bcount] == wall) {
            ins->mdotf[iface] = 0.0;
        }

        else if (ins->bcond[bcount] == inlet) {
            ins->mdotf[iface] = (ins->ubnd[bcount]*nx + ins->vbnd[bcount]*ny)*area;
        }

        else {

            dpdn = value[1]*nx + value[2]*ny;

            un = ins->uf[iface]*nx + ins->vf[iface]*ny;

            ins->mdotf[iface] = (un - ins->dt*dpdn)*area;
        }

        bcount++;
    }

    // 2) Correct velocity at cell centers of the mesh

    ins->max_div_vel = 0.0;

    for (icell = 0; icell < ins->mesh->n_cell; ++icell) {

        div_vel = 0.0; dp[0] = 0.0; dp[1] = 0.0;
        vol = ins->mesh->cell[icell].vol;

        for (int iface_local = 0; iface_local < ins->mesh->cell[icell].nf; ++iface_local) {

            iface = ins->mesh->cell[icell].f[iface_local];
            nx     = ins->mesh->face[iface].n[0];
            ny     = ins->mesh->face[iface].n[1];
            area   = ins->mesh->face[iface].area;
            snsign = ins->mesh->cell[icell].sn[iface_local];

            dp[0] += ins->pf[iface]*nx*area*snsign;
            dp[1] += ins->pf[iface]*ny*area*snsign;
            div_vel += ins->mdotf[iface]*snsign;
        }

        div_vel /= vol;

        dp[0] /= vol; dp[1] /= vol;

        ins->u[icell] -= ins->dt*dp[0]; ins->v[icell] -= ins->dt*dp[1];
        ins->divU[icell] = div_vel;

        if (fabs(div_vel) > ins->max_div_vel)
            ins->max_div_vel = fabs(div_vel);
    }
}

//---------------------------------------------------------------------
// Plot the solution (cell data) in ASCII VTK format
//---------------------------------------------------------------------

void ins_solver_plot_vtk(ins_solver* ins) {

    int icell;
    char filename[20];

    sprintf(filename, "sol-%08d.vtk", ins->timestep);

    triangulation_plot_vtk(ins->mesh, filename); // Plot mesh connectivity

    FILE *file_ptr;

    file_ptr = fopen(filename, "a");

    if (file_ptr == NULL)
        printf("Error opening the file.\n");

    fprintf(file_ptr, "CELL_DATA %d\n", ins->mesh->n_cell);

    // Plot pressure

    fprintf(file_ptr, "SCALARS Pressure double 1\n");
    fprintf(file_ptr, "LOOKUP_TABLE default\n");

    for (icell = 0; icell < ins->mesh->n_cell; ++icell)
        fprintf(file_ptr, "%.5e\n", ins->p[icell]);

    fprintf(file_ptr, "SCALARS divU double 1\n");
    fprintf(file_ptr, "LOOKUP_TABLE default\n");

    for (icell = 0; icell < ins->mesh->n_cell; ++icell)
        fprintf(file_ptr, "%.5e\n", ins->divU[icell]);


    // velocity

    fprintf(file_ptr, "VECTORS Velocity double\n");

    for (icell = 0; icell < ins->mesh->n_cell; ++icell)
        fprintf(file_ptr, "%.5e %.5e %.5e\n", ins->u[icell], ins->v[icell], 0.0);

    fclose(file_ptr);
}

//---------------------------------------------------------------------
// Calculate residue for monitoring convergence
//---------------------------------------------------------------------

void ins_solver_calc_residue(ins_solver* ins) {

    int icell;
    double V, V0;

    ins->residual = 0.0;

    for (icell = 0; icell < ins->mesh->n_cell; ++icell) {
        V0 = hypot(ins->u0[icell], ins->v0[icell]);
        V  = hypot(ins->u[icell], ins->v[icell]);

        ins->residual += (V-V0)*(V-V0);
    }

    ins->residual = sqrt(ins->residual)/(ins->dt*(double)ins->mesh->n_cell);
}

//---------------------------------------------------------------------
// Compute lift and drag coefficients for boundary with tag btag
//---------------------------------------------------------------------

void ins_solver_calc_force_coeffs(ins_solver* ins, int btag, double uref, double lref,  double* cl, double* cd) {

    int iface, bcount;
    double dudx, dudy, dvdx, dvdy;
    double nx, ny, area, value[ndof], p;
    double fl, fd;

    fl = 0.0; fd = 0.0;

    // 1b) -> Boundary faces

    bcount = 0;

    for (iface = ins->mesh->n_int_face; iface < ins->mesh->n_face; ++iface) {

        if (ins->mesh->face[iface].btag == btag) {

            area = ins->mesh->face[iface].area;
            nx   = ins->mesh->face[iface].n[0]; ny = ins->mesh->face[iface].n[1];

            evaluate_boundary_face(iface, ins->mesh, &ins->clsq_bnd_vel[bcount], ins->u, ins->ubnd[bcount], value);
            dudx = value[1]; dudy = value[2];

            evaluate_boundary_face(iface, ins->mesh, &ins->clsq_bnd_vel[bcount], ins->v, ins->vbnd[bcount], value);
            dvdx = value[1]; dvdy = value[2];

            evaluate_boundary_face(iface, ins->mesh, &ins->clsq_bnd_prs[bcount], ins->p, ins->pbnd[bcount], value);
            p = value[0];

            fd += -(-p*nx + ins->nu*(2.*nx*dudx + ny*(dvdx+dudy)))*area;
            fl += -(-p*ny + ins->nu*(2.*ny*dvdy + nx*(dvdx+dudy)))*area;
        }

        bcount++;
    }

    *cd = 2.*fd/(uref*uref*lref);
    *cl = 2.*fl/(uref*uref*lref);
}

//---------------------------------------------------------------------
// Put everything together and run the problem
//---------------------------------------------------------------------

void ins_solver_run(ins_solver* ins) {

    ins_solver_compute_reconstruction_data(ins); // Compute reconstruction data
    ins_solver_assemble_ppe_matrix(ins);         // Assemble matrix for solving pressure poisson equation

    ins_solver_plot_vtk(ins);

    while (ins->time < ins->tend) {

        ins_solver_calc_time_step(ins);

        ins_solver_calc_auxiliary_velocity(ins);

        ins_solver_calc_div_vel(ins);

        ins_solver_solve_ppe(ins);

        ins_solver_correct_velocity(ins);

        ins->time += ins->dt;
        ins->timestep++;

        if (ins->timestep%5 == 0) {
            ins_solver_calc_residue(ins);
            printf("step = %d, time = %f\n", ins->timestep, ins->time);
            printf("no. of PPE iterations = %d\n", ins->n_ppe_iter);
            printf("max-div-vel = %.3e, residue = %.3e\n", ins->max_div_vel, ins->residual);
            printf("------------------------------------------------------\n");

            if (ins->residual < 1.0e-6) break;
        }

        if (ins->timestep%500 == 0)
            ins_solver_plot_vtk(ins);
    }

    ins_solver_plot_vtk(ins);
}



