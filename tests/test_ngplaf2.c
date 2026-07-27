#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <lafp/gmlib.h>
#include <lafp/matrix_utils.h>
#include <lafp/mdarray.h>
#include <lafp/NGPlaf2.h>

/**
 * @file test_ngplaf2.c
 * @brief Integration test for NGPlaf2 multivariate LAFP sampler (Durante et al. 2014).
 */

int main(void) {
    int Nt = 50;
    int Np = 3;
    int NK = 2;
    int NL = 2;
    int Niter = 50;
    int Nsamps = Niter / 5;

    // Allocate synthetic observation matrix y (50 x 3) and tobs vector (50)
    gsl_matrix *y = gsl_matrix_alloc(Nt, Np);
    gsl_vector *tobs = gsl_vector_alloc(Nt);

    for (int i = 0; i < Nt; i++) {
        double t_val = (double)i * 0.1;
        gsl_vector_set(tobs, i, t_val);
        gsl_matrix_set(y, i, 0, sin(t_val));
        gsl_matrix_set(y, i, 1, cos(t_val));
        gsl_matrix_set(y, i, 2, sin(t_val) + 0.5 * cos(t_val));
    }

    // Allocate prior vectors
    gsl_vector *sigPrior = gsl_vector_alloc(2);
    gsl_vector *epsPrior = gsl_vector_alloc(2);
    gsl_vector *ksiPrior = gsl_vector_alloc(2);
    gsl_vector *APrior   = gsl_vector_alloc(2);
    gsl_vector *psiPrior = gsl_vector_alloc(2);
    gsl_vector *BPrior   = gsl_vector_alloc(2);
    gsl_vector *aaPrior  = gsl_vector_alloc(2);

    gsl_vector_set(sigPrior, 0, 4.0);
    gsl_vector_set(sigPrior, 1, 4.0);
    gsl_vector_set(epsPrior, 0, 1.0);
    gsl_vector_set(epsPrior, 1, 0.1);
    gsl_vector_set(ksiPrior, 0, 1.0);
    gsl_vector_set(ksiPrior, 1, 1.0);
    gsl_vector_set(APrior,   0, 1.0);
    gsl_vector_set(APrior,   1, 1.0);
    gsl_vector_set(psiPrior, 0, 1.0);
    gsl_vector_set(psiPrior, 1, 1.0);
    gsl_vector_set(BPrior,   0, 1.0);
    gsl_vector_set(BPrior,   1, 1.0);
    gsl_vector_set(aaPrior,  0, 3.0);
    gsl_vector_set(aaPrior,  1, 3.0);

    // Allocate output containers
    marray3d *Ksi_out = marray3d_alloc(Nsamps, Nt, NK * NL);
    marray3d *Psi_out = marray3d_alloc(Nsamps, Nt, NK);
    marray3d *yhat_out = marray3d_alloc(Nsamps, Nt, Np);
    marray3d *mu_out = marray3d_alloc(Nsamps, Nt, Np);
    marray4d *Sigma_out = marray4d_alloc(Nsamps, Nt, Np, Np);

    NGPlaf2 *laf = NGPlaf2_New();
    if (laf == NULL) GMERR(-10);

    if (NGPlaf2_construct(laf, tobs, y, NK, NL, Niter,
                          Ksi_out, Psi_out, yhat_out, mu_out, Sigma_out,
                          sigPrior, epsPrior, ksiPrior, APrior,
                          psiPrior, BPrior, aaPrior)) {
        GMERR(-20);
    }

    if (NGPlaf2_operations(laf)) GMERR(-30);

    // Verify non-NaN output values
    double mu_val = marray3d_get(mu_out, 0, 0, 0);
    double sig_val = marray4d_get(Sigma_out, 0, 0, 0, 0);

    if (mu_val != mu_val || sig_val != sig_val) {
        fprintf(stderr, "Error: NGPlaf2 produced NaN values.\n");
        GMERR(-40);
    }

    printf("[TEST_NGPLAF2] Successfully executed 50 iterations, Nsamps=%d.\n", Nsamps);
    printf("  mu_out[0,0,0]    = %.6f\n", mu_val);
    printf("  Sigma_out[0,0,0,0] = %.6f\n", sig_val);

    // Resource cleanup
    NGPlaf2_free(laf);
    gsl_matrix_free(y);
    gsl_vector_free(tobs);
    gsl_vector_free(sigPrior);
    gsl_vector_free(epsPrior);
    gsl_vector_free(ksiPrior);
    gsl_vector_free(APrior);
    gsl_vector_free(psiPrior);
    gsl_vector_free(BPrior);
    gsl_vector_free(aaPrior);
    marray3d_free(Ksi_out);
    marray3d_free(Psi_out);
    marray3d_free(yhat_out);
    marray3d_free(mu_out);
    marray4d_free(Sigma_out);

    return 0;

GMERRH("main", 1);
}
