#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <lafp/gmlib.h>
#include <lafp/mdarray.h>
#include <lafp/KalmanFilter2.h>

/**
 * @file test_kalman.c
 * @brief Unit test for Kalman2 state-space forward filtering.
 */

int main(void) {
    int Nt = 2;
    int Np = 1;
    int Nm = 1;
    int Nr = 1;

    // Allocate inputs
    gsl_matrix *y = gsl_matrix_alloc(Nt, Np);
    gsl_vector *a_init = gsl_vector_alloc(Nm);
    gsl_matrix *P_init = gsl_matrix_alloc(Nm, Nm);

    marray3d *Z = marray3d_alloc(Nt, Np, Nm);
    marray3d *H = marray3d_alloc(Nt, Np, Np);
    marray3d *T = marray3d_alloc(Nt, Nm, Nm);
    marray3d *R = marray3d_alloc(Nt, Nm, Nr);
    marray3d *Q = marray3d_alloc(Nt, Nr, Nr);

    // Allocate outputs
    gsl_matrix *v = gsl_matrix_alloc(Nt, Np);
    marray3d *K = marray3d_alloc(Nt, Nm, Np);
    marray3d *Finv = marray3d_alloc(Nt, Np, Np);

    // Populate model parameters (1D random walk + noise)
    gsl_matrix_set(y, 0, 0, 1.2);
    gsl_matrix_set(y, 1, 0, 1.8);

    gsl_vector_set(a_init, 0, 0.0);
    gsl_matrix_set(P_init, 0, 0, 1.0);

    marray3d_set_constant(Z, 1.0);
    marray3d_set_constant(H, 0.5);
    marray3d_set_constant(T, 1.0);
    marray3d_set_constant(R, 1.0);
    marray3d_set_constant(Q, 0.2);

    // Instantiate and configure Kalman filter
    Kalman2 *kFilter = Kalman2_New();
    if (kFilter == NULL) GMERR(-10);

    if (Kalman2_construct(kFilter, y, a_init, P_init, Z, H, T, R, Q, v, K, Finv)) {
        GMERR(-20);
    }

    // Run forward filter operations
    if (Kalman2_operations(kFilter)) GMERR(-30);

    // Basic assertion checks on residuals
    double v0 = gsl_matrix_get(v, 0, 0);
    if (v0 <= 0.0) GMERR(-40);

    // Resource cleanup
    Kalman2_free(kFilter);
    gsl_matrix_free(y);
    gsl_vector_free(a_init);
    gsl_matrix_free(P_init);
    marray3d_free(Z);
    marray3d_free(H);
    marray3d_free(T);
    marray3d_free(R);
    marray3d_free(Q);
    gsl_matrix_free(v);
    marray3d_free(K);
    marray3d_free(Finv);

    return 0;

GMERRH("main", 1);
}
