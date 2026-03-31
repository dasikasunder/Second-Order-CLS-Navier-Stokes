/*
 * triangulation.c
 *      Author: sunder
 */

#include "triangulation.h"

// Additional helper functions

void form_vrtx_data(triangulation*);
void form_face_data(triangulation*, int*, int*, int*);
void form_cell_data(triangulation*);
void calc_cell_metrics(triangulation*);
void calc_face_metrics(triangulation*);

double triangle_area(double x1, double x2, double x3, double y1, double y2, double y3) {
    return 0.5*( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );
}

// Read mesh from the mesh file

triangulation* triangulation_allocate(const char* filename) {

    int i, j, rval, tmp, counter;
    int *bvrtx1, *bvrtx2, *btag, tag, n_tag;
    char line[50];

    triangulation* tria = (triangulation*)malloc(sizeof(triangulation));

    /* 1) Open and read the mesh file */

    FILE *file_ptr;
    file_ptr = fopen(filename, "r");

    if (file_ptr == NULL) {
        printf("Error: Could not open file %s.\n", filename);
        exit(EXIT_FAILURE);
    }

    // 1a) Get dimensions of the mesh. Must be two

    rval = fscanf(file_ptr, "%s %d", line, &tria->n_dim);

    if (tria->n_dim != 2) {
        printf("No. of dimensions in mesh is %d.\n", tria->n_dim);
        printf("Only 2D meshes are supported\n");
        exit(EXIT_FAILURE);
    }

    // 1b) Get the number of cells in the mesh and allocate memory for cell array

    rval = fscanf(file_ptr, "%s %d", line, &tria->n_cell);
    tria->cell = (Cell *) malloc(tria->n_cell*sizeof(Cell));

    if (tria->cell == NULL) {
        printf("Memory allocation to cell array failed!\n");
        exit(EXIT_FAILURE);
    }

    tria->n_tria = 0; tria->n_quad = 0;

    // 1c) Read mesh connectivity

    for (i = 0; i < tria->n_cell; ++i) {

        rval = fscanf(file_ptr, "%d\n", &counter);

        if (counter == 5) {      // Triangle

            rval = fscanf(file_ptr, "%d %d %d %d", &tria->cell[i].v[0],
                          &tria->cell[i].v[1], &tria->cell[i].v[2], &tmp);

            tria->cell[i].nv = 3;
            tria->cell[i].nf = 3;
            tria->n_tria++;

        }

        else if (counter == 9) { // Quadrilateral

            rval = fscanf(file_ptr, "%d %d %d %d %d", &tria->cell[i].v[0],
                          &tria->cell[i].v[1], &tria->cell[i].v[2],
                          &tria->cell[i].v[3], &tmp);

            tria->cell[i].nv = 4; // This value will be filled later
            tria->cell[i].nf = 4;
            tria->n_quad++;
        }

        else {
            printf("Element type %d is not supported.\n", tria->cell[i].nv);
            exit(EXIT_FAILURE);
        }
    }

    // 1d) Get number of vertices in the mesh

    rval = fscanf(file_ptr, "%s %d", line, &tria->n_vrtx);

    tria->vrtx = (Vrtx *) malloc(tria->n_vrtx*sizeof(Vrtx));

    if (tria->vrtx == NULL) {
        printf("Memory allocation to vrtx array failed!\n");
        exit(EXIT_FAILURE);
    }

    // 1d) Read coordinates of the mesh

    for (i = 0; i < tria->n_vrtx; ++i)
        rval = fscanf(file_ptr, "%lf %lf %d\n", &tria->vrtx[i].x[0], &tria->vrtx[i].x[1], &tmp);

    // 1e) Read boundary condition data

    rval = fscanf(file_ptr, "%s %d", line, &tria->n_bnd_tag);

    tria->n_bnd_face = 0;

    bvrtx1 = (int *)malloc(tria->n_vrtx * sizeof(int));
    bvrtx2 = (int *)malloc(tria->n_vrtx * sizeof(int));
    btag   = (int *)malloc(tria->n_vrtx * sizeof(int));

    counter = 0;

    for (i = 0; i < tria->n_bnd_tag; ++i) {

        rval = fscanf(file_ptr, "%s %d", line, &tag);
        rval = fscanf(file_ptr, "%s %d", line, &n_tag);

        tria->n_bnd_face += n_tag;

        for (j = 0; j < n_tag; ++j) {

            rval = fscanf(file_ptr, "%d %d %d\n", &tmp, &bvrtx1[counter], &bvrtx2[counter]);

            btag[counter] = tag;
            counter++;
        }
    }

    (void)rval;
    fclose(file_ptr);

    /* 2) Form additional connectivity */

    form_vrtx_data(tria);                          // 2a) Vertex to cell
    form_face_data(tria, bvrtx1, bvrtx2, btag);    // 2b) Find face neighbours of a cell
    form_cell_data(tria);                          // 2c) Form cell vertex and cell face data

    free(bvrtx1); free(bvrtx2); free(btag);

    /* 3) Calculate grid metrics */

    calc_cell_metrics(tria); // 3a) Cell metrics
    calc_face_metrics(tria); // 3b) Face metrics

    return tria;
}

/* Form vertex data */

void form_vrtx_data(triangulation* tria) {

    int i, k, vk;

    for (i = 0; i < tria->n_vrtx; ++i) // Initialize to zero
        tria->vrtx[i].nc = 0;

    // Count the # of cells around a vertex

    for (i = 0; i < tria->n_cell; ++i) {
        for (k = 0; k < tria->cell[i].nv; ++k) {
            vk = tria->cell[i].v[k];
            tria->vrtx[vk].nc = tria->vrtx[vk].nc + 1;
        }
    }

    // Allocate memory for vertex to cell array

    for (i = 0; i < tria->n_vrtx; ++i)
        allocate_1d(tria->vrtx[i].c, tria->vrtx[i].nc);

    // Re-initialize the # of cells around a vertex

    for (i = 0; i < tria->n_vrtx; ++i)
        tria->vrtx[i].nc = 0;

    for (i = 0; i < tria->n_cell; ++i) {
        for (k = 0; k < tria->cell[i].nv; ++k) {
            vk = tria->cell[i].v[k];
            tria->vrtx[vk].c[tria->vrtx[vk].nc] = i;
            tria->vrtx[vk].nc = tria->vrtx[vk].nc + 1;
        }
    }
}

/* Find face neighbours of a cell */

void form_face_data(triangulation* tria, int* bvrtx1, int* bvrtx2, int* btag) {

    int i, j, k, kk, v1, v2, ck, vk1, vk2, vk3, vk4;
    int ***commonv;

    allocate_3d(commonv,tria->n_cell, 4, 2);

    for (i = 0; i < tria->n_cell; ++i) { // Loop over all the cells

        tria->cell[i].nn = 0; // Initialize the # of neighbors around a cell i.

        for (j = 0; j < 4; ++j)
            tria->cell[i].n[j] = -1; // Initialize neighbor list to -1

        // Loop over the vertices of the cell

        for (k = 0; k < tria->cell[i].nv; ++k) {

            // Define a face by two vertices [v1,v2], for which a neighbor is sought.

            v1 = tria->cell[i].v[k];

            if (k+1 >= tria->cell[i].nv)
                v2 = tria->cell[i].v[0];
            else
                v2 = tria->cell[i].v[k+1];

            // Look for a neighbor in the cells around v1

            for (kk = 0; kk < tria->vrtx[v1].nc; ++kk) {

                ck = tria->vrtx[v1].c[kk];

                if (tria->cell[ck].nv == 3) { // Neighbor cell is a triangle

                    vk1 = tria->cell[ck].v[0];
                    vk2 = tria->cell[ck].v[1];
                    vk3 = tria->cell[ck].v[2];

                    // Check if any of the face matches (v1,v2)

                    if ( (v1 == vk1 && v2 == vk2) ||
                         (v1 == vk2 && v2 == vk1) ||
                         (v1 == vk2 && v2 == vk3) ||
                         (v1 == vk3 && v2 == vk2) ||
                         (v1 == vk1 && v2 == vk3) ||
                         (v1 == vk3 && v2 == vk1) ) {

                        if (ck != i) {
                            tria->cell[i].n[tria->cell[i].nn] = ck;
                            commonv[i][tria->cell[i].nn][0] = v1;
                            commonv[i][tria->cell[i].nn][1] = v2;
                            tria->cell[i].nn++;
                        }
                    }
                }

                else { // Neighbor cell is a quadrilateral

                    vk1 = tria->cell[ck].v[0];
                    vk2 = tria->cell[ck].v[1];
                    vk3 = tria->cell[ck].v[2];
                    vk4 = tria->cell[ck].v[3];

                    // Check if any of the face matches (v1,v2)

                    if ( (v1 == vk1 && v2 == vk2) ||
                         (v1 == vk2 && v2 == vk1) ||
                         (v1 == vk2 && v2 == vk3) ||
                         (v1 == vk3 && v2 == vk2) ||
                         (v1 == vk4 && v2 == vk3) ||
                         (v1 == vk3 && v2 == vk4) ||
                         (v1 == vk1 && v2 == vk4) ||
                         (v1 == vk4 && v2 == vk1) ) {


                        if (ck != i) {
                            tria->cell[i].n[tria->cell[i].nn] = ck;
                            commonv[i][tria->cell[i].nn][0] = v1;
                            commonv[i][tria->cell[i].nn][1] = v2;
                            tria->cell[i].nn++;
                        }
                    }
                }
            } // vertex neighbour loop
        } // Vertex loop
    } // Cell loop

    // First, count the number of total interior faces

    tria->n_face = 0;

    for (i = 0; i < tria->n_cell; ++i) {
        for (k = 0; k < tria->cell[i].nn; ++k) {
            // This identifies a face and avoids double count.
            if ( i < tria->cell[i].n[k] )
                tria->n_face++;
        }
    }

    tria->n_int_face = tria->n_face;


    tria->n_face = tria->n_face+tria->n_bnd_face;

    // Allocate and fill interior face arrays

    allocate_1d(tria->face, tria->n_face);

    // We need to re-count 'n_face'. So, initialize it.

    tria->n_face = 0;

    for (i = 0; i < tria->n_cell; i++) {
        for (k = 0; k < tria->cell[i].nn; ++k) {
            // This identifies a face and avoids double count.
            if ( i < tria->cell[i].n[k] ) {

                tria->face[tria->n_face].c[0] = i;
                tria->face[tria->n_face].c[1] = tria->cell[i].n[k];
                tria->face[tria->n_face].v[0] = commonv[i][k][0];
                tria->face[tria->n_face].v[1] = commonv[i][k][1];
                tria->face[tria->n_face].btag = 0;
                tria->face[tria->n_face].at_boundary = false;
                tria->n_face++;
            }
        }
    }

    // Add boundary faces

    for (i  = 0; i < tria->n_bnd_face; ++i) {

        tria->face[tria->n_face].v[0] = bvrtx1[i];
        tria->face[tria->n_face].v[1] = bvrtx2[i];
        tria->face[tria->n_face].btag = btag[i];
        tria->face[tria->n_face].at_boundary = true;
        tria->n_face++;
    }

    // Add neighbouring cell to boundary faces

    for (i = tria->n_face - tria->n_bnd_face; i < tria->n_face; ++i) {

        v1 = tria->face[i].v[0];
        v2 = tria->face[i].v[1];

        for (kk = 0; kk < tria->vrtx[v1].nc; ++kk) {

            ck = tria->vrtx[v1].c[kk];

            if (tria->cell[ck].nv == 3) { // Neighbor cell is a triangle

                vk1 = tria->cell[ck].v[0];
                vk2 = tria->cell[ck].v[1];
                vk3 = tria->cell[ck].v[2];

                // Check if any of the face matches (v1,v2)

                if ( (v1 == vk1 && v2 == vk2) ||
                    (v1 == vk2 && v2 == vk1) ||
                    (v1 == vk2 && v2 == vk3) ||
                    (v1 == vk3 && v2 == vk2) ||
                    (v1 == vk1 && v2 == vk3) ||
                    (v1 == vk3 && v2 == vk1) ) {

                    tria->face[i].c[0] = ck;
                    tria->face[i].c[1] = -1;
                    break;

                }
            }

            else { // Neighbor cell is a quadrilateral

                vk1 = tria->cell[ck].v[0];
                vk2 = tria->cell[ck].v[1];
                vk3 = tria->cell[ck].v[2];
                vk4 = tria->cell[ck].v[3];

                // Check if any of the face matches (v1,v2)

                if ( (v1 == vk1 && v2 == vk2) ||
                     (v1 == vk2 && v2 == vk1) ||
                     (v1 == vk2 && v2 == vk3) ||
                     (v1 == vk3 && v2 == vk2) ||
                     (v1 == vk4 && v2 == vk3) ||
                     (v1 == vk3 && v2 == vk4) ||
                     (v1 == vk1 && v2 == vk4) ||
                     (v1 == vk4 && v2 == vk1) ) {

                    tria->face[i].c[0] = ck;
                    tria->face[i].c[1] = -1;
                    break;
                }
            }
        } // vertex neighbour loop
    }

    // Free commonv array

    free_3d(commonv);
}

/* Form cell data */

void form_cell_data(triangulation* tria) {

    int i, j, k, ii, c1, c2, vk, kc, jcell, f, count, fstore[4];
    int candidate_cell, c[25], nc;
    bool already_added;

    // Create cell to face connectivity

    for (i = 0; i < tria->n_cell; ++i) // Initialize to zero
        tria->cell[i].nf = 0;

    // Populate cell->face list

    for (i = 0; i < tria->n_face; ++i) {
        c1 = tria->face[i].c[0];
        tria->cell[c1].f[tria->cell[c1].nf] = i;
        tria->cell[c1].nf++;

        if (! tria->face[i].at_boundary) {
            c2 = tria->face[i].c[1];
            tria->cell[c2].f[tria->cell[c2].nf] = i;
            tria->cell[c2].nf++;
        }
    }

    // Align faces and neighbour cell

    for (i = 0; i < tria->n_cell; ++i) {

        count = tria->cell[i].nf-1;

        for (j = 0; j < tria->cell[i].nf; ++j) {

            f = tria->cell[i].f[j];

            if (tria->face[f].c[0] == i)
                jcell = tria->face[f].c[1];
            else
                jcell = tria->face[f].c[0];

            for (k = 0; k < tria->cell[i].nf; ++k) {
                if (tria->cell[i].n[k] == jcell && jcell != -1) {
                    fstore[k] = f;
                }
            }

            if (jcell == -1) {
                fstore[count] = f;
                count--;
            }
        }

        for (j = 0; j < tria->cell[i].nf; ++j)
            tria->cell[i].f[j] = fstore[j];
    }

    // Find cells around a face

    for (i = 0; i < tria->n_face; ++i) {

        nc = 0;

        for (j = 0; j < 25; ++j)
            c[j] = -1;

        for (k = 0; k < 2; ++k) { // 2 is for two vertices of the face

            vk = tria->face[i].v[k];

            for (kc = 0; kc < tria->vrtx[vk].nc; ++kc) {

                candidate_cell = tria->vrtx[vk].c[kc];

                already_added = false;

                if (nc > 0) {
                    for (ii = 0; ii < nc; ++ii) {
                        if ( candidate_cell == c[ii] )
                            already_added = true;
                    }
                }

                if (already_added == false &&
                    candidate_cell != tria->face[i].c[0] &&
                    candidate_cell != tria->face[i].c[1]) {
                    c[nc] = candidate_cell;
                    nc++;
                }
            }
        }

        tria->face[i].nvnb = nc;
        allocate_1d(tria->face[i].vnb, nc);
        for (k = 0; k < nc; ++k)
            tria->face[i].vnb[k] = c[k];

    } // Face loop
}

/* Calculate centroids, volumes and other metrics related to cells */

void calc_cell_metrics(triangulation* tria) {

    int i, v1, v2, v3, v4;
    double x1, x2, x3, x4, y1, y2, y3, y4;


    for (i = 0; i < tria->n_cell; ++i) {

        if (tria->cell[i].nf == 3) { // Triangle

            v1 = tria->cell[i].v[0]; v2 = tria->cell[i].v[1];
            v3 = tria->cell[i].v[2];

            x1 = tria->vrtx[v1].x[0]; y1 = tria->vrtx[v1].x[1];
            x2 = tria->vrtx[v2].x[0]; y2 = tria->vrtx[v2].x[1];
            x3 = tria->vrtx[v3].x[0]; y3 = tria->vrtx[v3].x[1];

            tria->cell[i].vol = triangle_area(x1,x2,x3,y1,y2,y3);
            tria->cell[i].x[0] = (1./3.)*(x1+x2+x3);
            tria->cell[i].x[1] = (1./3.)*(y1+y2+y3);
        }

        else { // Quadrilateral

            v1 = tria->cell[i].v[0]; v2 = tria->cell[i].v[1];
            v3 = tria->cell[i].v[2]; v4 = tria->cell[i].v[3];

            x1 = tria->vrtx[v1].x[0]; y1 = tria->vrtx[v1].x[1];
            x2 = tria->vrtx[v2].x[0]; y2 = tria->vrtx[v2].x[1];
            x3 = tria->vrtx[v3].x[0]; y3 = tria->vrtx[v3].x[1];
            x4 = tria->vrtx[v4].x[0]; y4 = tria->vrtx[v4].x[1];

            tria->cell[i].vol = triangle_area(x1,x2,x4,y1,y2,y4) + triangle_area(x2,x3,x4,y2,y3,y4);
            tria->cell[i].x[0] = (1./4.)*(x1+x2+x3+x4);
            tria->cell[i].x[1] = (1./4.)*(y1+y2+y3+y4);
        }
    }
}

/* Calculate centers, areas and other metrics related to faces */

void calc_face_metrics(triangulation* tria) {

	int i, j, k, v1, v2, c1, c2;
    double x1, x2, y1, y2;
    double xc1, yc1, xc2, yc2;

    for (i = 0; i < tria->n_face; ++i) {

    	v1 = tria->face[i].v[0]; v2 = tria->face[i].v[1];
        x1 = tria->vrtx[v1].x[0]; y1 = tria->vrtx[v1].x[1];
        x2 = tria->vrtx[v2].x[0]; y2 = tria->vrtx[v2].x[1];

        tria->face[i].x[0] = 0.5*(x1+x2); tria->face[i].x[1] = 0.5*(y1+y2);
        tria->face[i].area = hypot(x2-x1,y2-y1);
        tria->face[i].n[0] = -(y2-y1)/tria->face[i].area;
        tria->face[i].n[1] =  (x2-x1)/tria->face[i].area;

        c1 = tria->face[i].c[0];
        xc1 = tria->cell[c1].x[0]; yc1 = tria->cell[c1].x[1];

        if (tria->face[i].at_boundary) {
            xc2 = tria->face[i].x[0]; yc2 = tria->face[i].x[1];
        }
        else {
            c2 = tria->face[i].c[1];
            xc2 = tria->cell[c2].x[0]; yc2 = tria->cell[c2].x[1];
        }

        tria->face[i].r = hypot(xc2-xc1,yc2-yc1);
        tria->face[i].e[0] = (xc2-xc1)/tria->face[i].r;
        tria->face[i].e[1] = (yc2-yc1)/tria->face[i].r;

        // Make sure n and e are aligned in the same direction

        if (tria->face[i].e[0]*tria->face[i].n[0] + tria->face[i].e[1]*tria->face[i].n[1] < 0.0) {
            tria->face[i].n[0] = -tria->face[i].n[0];
            tria->face[i].n[1] = -tria->face[i].n[1];
        }
    }

    // Add proper normal sign to cell faces

    for (i = 0; i < tria->n_cell; ++i) {

    	for (j = 0; j < tria->cell[i].nf; ++j) {

    		k = tria->cell[i].f[j];

    		if (tria->face[k].c[0] == i)
    			tria->cell[i].sn[j] = 1.0;
    		else
    			tria->cell[i].sn[j] = -1.0;
    	}
    }
}

/* Verify triangulation data */

void triangulation_verify(const triangulation* tria) {

	int icell, iface, result;
	const double rel_accuracy = 1.0e-10;
	double vol1 = 0.0, vol2 = 0.0;

	// 1) Total volume of the triangulation

	// 1a) Total volume calculated by adding volumes of the cells

	for (icell = 0; icell < tria->n_cell; ++icell)
		vol1 += tria->cell[icell].vol;

	// 1b) Total volume calculated using Green-Gauss theorem for boundary faces

	for (iface = tria->n_int_face; iface < tria->n_face; ++iface)

		vol2 += 0.5*(tria->face[iface].x[0]*tria->face[iface].n[0] +
				     tria->face[iface].x[1]*tria->face[iface].n[1])*tria->face[iface].area;

	result = gsl_fcmp(vol1, vol2, rel_accuracy);

	if (result == 0)
		printf("Domain Volume Test: PASS\n");
	else
		printf("Domain Volume Test: FAIL\n");

	// 2) Boundary normal test (sum should be zero)

}

/* Plot triangulation in VTK format */

void triangulation_plot_vtk(const triangulation* tria, const char* filename) {

    int i, j, cell_size = 0;

    FILE *fptr;
    fptr = fopen(filename, "w");

    // Write header

    fprintf(fptr,"# vtk DataFile Version 2.0\n");
    fprintf(fptr,"Generic Title\n");
    fprintf(fptr,"ASCII\n");
    fprintf(fptr,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(fptr, "POINTS %d double\n", tria->n_vrtx);

    // 1) Write vertex info

    for (i = 0; i < tria->n_vrtx; ++i)
        fprintf(fptr, "%f %f %f\n", tria->vrtx[i].x[0], tria->vrtx[i].x[1], 0.0);

    // 2) Write cell connectivity info

    for (i = 0; i < tria->n_cell; ++i)
        cell_size += tria->cell[i].nv + 1;

    fprintf(fptr, "CELLS %d %d\n", tria->n_cell, cell_size);

    for (i = 0; i < tria->n_cell; ++i) {
        fprintf(fptr, "%d ", tria->cell[i].nv);
        for (j = 0; j < tria->cell[i].nv; ++j)
            fprintf(fptr, "%d ", tria->cell[i].v[j]);
        fprintf(fptr,"\n");
    }

    fprintf(fptr, "CELL_TYPES %d\n", tria->n_cell);

    for (i = 0; i < tria->n_cell; ++i) {
        if (tria->cell[i].nv == 3) // Triangle
            fprintf(fptr, "%d\n",5);
        else // Quadrilateral
            fprintf(fptr, "%d\n",9);
    }

     fclose(fptr);
}

/* Release all the memory allocated to triangulation object */

void triangulation_free(triangulation* tria) {

    int i;

    for (i = 0; i < tria->n_vrtx; ++i)
        free(tria->vrtx[i].c);

    for (i = 0; i < tria->n_face; ++i)
        free(tria->face[i].vnb);

    free(tria->vrtx);
    free(tria->cell);
    free(tria->face);

    free(tria);
}


