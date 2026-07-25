#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <lafp/gmlib.h>
#include <lafp/matrix_utils.h>
#include <lafp/mdarray.h>
#include <lafp/NGPmcmc.h>

/**
 * @file test_ngpmcmc.c
 * @brief Integration test for the Non-Gaussian Process (NGP) MCMC sampler.
 *
 * Reads input observation data (y.txt and tobs.txt), initializes model
 * hyperparameters, executes Niter MCMC sampling iterations, writes the
 * posterior output files, and verifies clean memory deallocation.
 */

int main(void) {
    FILE *fp;
    int Niter = 2000;
    int Nt = 1001;

    // Allocate observation matrix, hyperparameters, and posterior containers
    gsl_matrix *y = gsl_matrix_alloc(Nt, 1);
    gsl_vector *sigU = gsl_vector_alloc(1);
    gsl_vector *sigA = gsl_vector_alloc(1);
    gsl_vector *tobs = gsl_vector_alloc(Nt);
    marray3d *th = marray3d_alloc(Niter, Nt, 3);
    gsl_matrix *sig = gsl_matrix_alloc(Niter, 3);

    // Initial hyperparameter values
    double a = 1.0;
    double b = 1.0;
    double sigMu = 4.0;
    double sigAlph = 4.0;
    double sigEps = 2.0;

    gsl_vector_set(sigU, 0, 1000.0);
    gsl_vector_set(sigA, 0, 5.0);

    // Load observation values
    fp = fopen("y.txt", "r");
    if (fp == NULL) GMERR(-11);
    gsl_matrix_fscanf(fp, y);
    fclose(fp);

    // Load time observations
    fp = fopen("tobs.txt", "r");
    if (fp == NULL) GMERR(-21);
    gsl_vector_fscanf(fp, tobs);
    fclose(fp);

    // Construct and configure MCMC sampler object
    NGPmcmc *myMCMC = NGPmcmc_New();
    if (myMCMC == NULL) GMERR(-30);
    if (NGPmcmc_construct(myMCMC, y, tobs, Niter, sigU, sigA,
                          sigEps, sigMu, sigAlph, a, b, th, sig)) {
        GMERR(-31);
    }

    // Execute MCMC sampling loop
    if (NGPmcmc_operations(myMCMC)) GMERR(-41);

    // Write posterior sample output files (AMY_M1_Theta.txt, AMY_M1_Sig.txt)
    if (NGPmcmc_writeOutputs(myMCMC, "AMY_M1")) GMERR(-101);

    // Resource cleanup
    NGPmcmc_free(myMCMC);
    gsl_matrix_free(y);
    gsl_vector_free(sigU);
    gsl_vector_free(sigA);
    marray3d_free(th);
    gsl_matrix_free(sig);
    gsl_vector_free(tobs);

    return 0;

GMERRH("main", 1);
}
