#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <lafp/gmlib.h>
#include <lafp/matrix_utils.h>

/**
 * @file test_matrix_utils.c
 * @brief Unit tests for linear algebra matrix and vector utility functions.
 */

int main(void) {
    int N = 3;
    gsl_matrix *A = gsl_matrix_alloc(N, N);
    gsl_matrix *B = gsl_matrix_alloc(N, N);
    gsl_matrix *C = gsl_matrix_alloc(N, N);
    gsl_vector *v = gsl_vector_alloc(N);

    // Initialize test matrix A as 3x3 identity
    gsl_matrix_set_identity(A);

    // Test symmetricMatrix & positiveDefinite on identity matrix
    if (symmetricMatrix(A) != 0) GMERR(-10);
    if (positiveDefinite(A) != 0) GMERR(-11);

    // Test matrix norms on 3x3 identity matrix (Frobenius norm should be sqrt(3))
    double normF = 0.0;
    if (mnormF(A, &normF)) GMERR(-20);
    if (fabs(normF - sqrt(3.0)) > 1e-6) GMERR(-21);

    double norm1 = 0.0;
    if (mnorm1(A, &norm1)) GMERR(-22);
    if (fabs(norm1 - 1.0) > 1e-6) GMERR(-23);

    double norm00 = 0.0;
    if (mnorm00(A, &norm00)) GMERR(-24);
    if (fabs(norm00 - 1.0) > 1e-6) GMERR(-25);

    // Test positiveVector and positiveMatrix
    gsl_vector_set_all(v, 2.5);
    if (positiveVector(v) != 0) GMERR(-30);

    gsl_matrix_set_all(B, 1.5);
    if (positiveMatrix(B) != 0) GMERR(-31);

    // Test dot_mm matrix multiplication (A * B = B)
    if (dot_mm(A, B, C)) GMERR(-40);
    if (gsl_matrix_get(C, 0, 0) != 1.5) GMERR(-41);

    // Test dot_amaT transformation (A * B * A^T = B)
    if (dot_amaT(A, B, C)) GMERR(-50);
    if (gsl_matrix_get(C, 0, 0) != 1.5) GMERR(-51);

    // Clean up
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(C);
    gsl_vector_free(v);

    return 0;

GMERRH("main", 1);
}
