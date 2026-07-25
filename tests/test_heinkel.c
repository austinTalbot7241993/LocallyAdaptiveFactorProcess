#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <lafp/gmlib.h>
#include <lafp/matrix_utils.h>

/**
 * @file test_heinkel.c
 * @brief Unit test for hankel matrix construction from vector inputs.
 */

int main(void) {
    int N = 3;
    int N2 = 5;
    int i;
    double temp;

    gsl_vector *c = gsl_vector_alloc(N);
    gsl_vector *r = gsl_vector_alloc(N2);
    gsl_matrix *out = gsl_matrix_alloc(N, N2);

    for (i = 0; i < N; i++) {
        temp = 2.1 * i;
        gsl_vector_set(c, i, temp);
    }
    gsl_vector_set(r, 0, 10.0);
    gsl_vector_set(r, 1, 11.0);
    gsl_vector_set(r, 2, 12.0);
    gsl_vector_set(r, 3, 14.0);
    gsl_vector_set(r, 4, 15.0);

    // Build Hankel matrix
    if (hankel(c, r, out)) GMERR(-11);

    // Verify print output
    printGSLMatrix(out);

    // Free resources
    gsl_vector_free(c);
    gsl_vector_free(r);
    gsl_matrix_free(out);

    return 0;

GMERRH("main", 1);
}
