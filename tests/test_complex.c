#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <lafp/gmlib.h>
#include <lafp/mdarray_complex.h>
#include <lafp/KalmanFilter2c.h>
#include <lafp/fastStateSmoother2c.h>

/**
 * @file test_complex.c
 * @brief Unit test for complex-valued multidimensional arrays, KalmanFilter2c, and FSS2c.
 */

int main(void) {
    int Nt = 2;
    int Np = 1;
    int Nm = 1;
    int Nr = 1;

    // Test 3D complex multidimensional array allocation, set, and get
    marray3d_complex *arr = marray3d_complex_alloc(Nt, Np, Nm);
    if (arr == NULL) GMERR(-10);

    gsl_complex c_val = gsl_complex_rect(3.5, -2.1);
    marray3d_complex_set(arr, 0, 0, 0, c_val);

    gsl_complex c_got = marray3d_complex_get(arr, 0, 0, 0);
    if (GSL_REAL(c_got) != 3.5 || GSL_IMAG(c_got) != -2.1) GMERR(-11);

    marray3d_complex_free(arr);

    // Test KalmanFilter2c allocation, construction, and operations
    gsl_matrix_complex *y = gsl_matrix_complex_alloc(Nt, Np);
    gsl_vector_complex *a_init = gsl_vector_complex_alloc(Nm);
    gsl_matrix_complex *P_init = gsl_matrix_complex_alloc(Nm, Nm);

    marray3d_complex *Z = marray3d_complex_alloc(Nt, Np, Nm);
    marray3d_complex *H = marray3d_complex_alloc(Nt, Np, Np);
    marray3d_complex *T = marray3d_complex_alloc(Nt, Nm, Nm);
    marray3d_complex *R = marray3d_complex_alloc(Nt, Nm, Nr);
    marray3d_complex *Q = marray3d_complex_alloc(Nt, Nr, Nr);

    gsl_matrix_complex *v = gsl_matrix_complex_alloc(Nt, Np);
    marray3d_complex *K = marray3d_complex_alloc(Nt, Nm, Np);
    marray3d_complex *Finv = marray3d_complex_alloc(Nt, Np, Np);
    gsl_matrix_complex *alpha = gsl_matrix_complex_alloc(Nt, Nm);

    gsl_complex one = gsl_complex_rect(1.0, 0.0);
    gsl_matrix_complex_set_all(y, one);
    gsl_vector_complex_set_all(a_init, gsl_complex_rect(0.0, 0.0));
    gsl_matrix_complex_set_all(P_init, one);

    marray3d_complex_set_constant(Z, one);
    marray3d_complex_set_constant(H, gsl_complex_rect(0.5, 0.0));
    marray3d_complex_set_constant(T, one);
    marray3d_complex_set_constant(R, one);
    marray3d_complex_set_constant(Q, gsl_complex_rect(0.2, 0.0));

    Kalman2c *kf_c = Kalman2c_New();
    if (kf_c == NULL) GMERR(-20);

    if (Kalman2c_construct(kf_c, y, a_init, P_init, Z, H, T, R, Q, v, K, Finv)) GMERR(-21);
    if (Kalman2c_operations(kf_c)) GMERR(-22);

    // Test Fast State Smoother Complex (FSS2c)
    FSS2c *fss_c = FSS2c_New();
    if (fss_c == NULL) GMERR(-30);

    if (FSS2c_construct(fss_c, v, K, Finv, a_init, P_init, Z, T, R, Q, alpha)) GMERR(-31);
    if (FSS2c_operations(fss_c)) GMERR(-32);

    // Resource Cleanup
    FSS2c_free(fss_c);
    Kalman2c_free(kf_c);

    gsl_matrix_complex_free(y);
    gsl_vector_complex_free(a_init);
    gsl_matrix_complex_free(P_init);
    marray3d_complex_free(Z);
    marray3d_complex_free(H);
    marray3d_complex_free(T);
    marray3d_complex_free(R);
    marray3d_complex_free(Q);
    gsl_matrix_complex_free(v);
    marray3d_complex_free(K);
    marray3d_complex_free(Finv);
    gsl_matrix_complex_free(alpha);

    return 0;

GMERRH("main", 1);
}
