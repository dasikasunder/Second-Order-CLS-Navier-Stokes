/*
 * includes.h
 *      Author: sunder
 */

#ifndef INCLUDES_H_
#define INCLUDES_H_

// Standard C libraries
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

// GSL libraries
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

// LIS library for solving sparse linear equations

#include <lis.h>

extern const int ndof; // No of degrees of freedom (3 for second order)

// Utility functions

// Allocating memory for arrays

void *malloc_or_exit(size_t, const char*, int);

#define xmalloc(nbytes) malloc_or_exit((nbytes), __FILE__, __LINE__)

#define allocate_1d(v,n) ((v) = xmalloc((n) * sizeof *(v)))

#define free_1d(v) do { free(v); v = NULL; } while (0)

#define allocate_2d(a, imax, jmax) do { \
    int i; \
    allocate_1d(a, (imax) + 1); \
    for (i = 0; i < (imax); i++) \
        allocate_1d((a)[i], (jmax)); \
    (a)[imax] = NULL; \
} while (0)

#define free_2d(a) do { \
    if (a != NULL) { \
        int i; \
        for (i = 0; (a)[i] != NULL; i++) \
            free_1d((a)[i]); \
        free_1d(a); \
        a = NULL; \
    } \
} while (0)

#define allocate_3d(a, imax, jmax, kmax) do { \
    int i, j; \
    allocate_1d(a, (imax) + 1); \
    for (i = 0; i < (imax); i++) { \
        allocate_1d((a)[i], (jmax)+1); \
        for (j = 0; j < (jmax); j++) { \
            allocate_1d((a)[i][j], (kmax)); \
        } \
        (a)[i][jmax] = NULL; \
    } \
    (a)[imax] = NULL; \
} while(0)

#define free_3d(a) do { \
if (a != NULL) { \
    int i, j; \
    for (i = 0; (a)[i] != NULL; i++) { \
        for (j = 0; (a)[i][j] != NULL; j++ ) { \
            free_1d((a)[i][j]); \
        } \
        free_1d((a)[i]); \
    } \
    free_1d(a); \
    a = NULL; \
} \
} while (0)

#endif /* INCLUDES_H_ */
