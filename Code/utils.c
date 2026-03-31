/*
 * utils.c
 *      Author: sunder
 */

#include "includes.h"


const int ndof = 3;

/* Allocate memory or exit in case of failure */

void *malloc_or_exit(size_t nbytes, const char *file, int line) {
    void *x;

    if ((x = malloc(nbytes)) == NULL) {
        fprintf(stderr, "%s:line %d: malloc() of %zu bytes failed\n", file, line, nbytes);
        exit(EXIT_FAILURE);
    }

    else
        return x;
}
