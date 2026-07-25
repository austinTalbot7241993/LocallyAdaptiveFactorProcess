#ifndef _ABT_NGPMCMC_H
#define _ABT_NGPMCMC_H

/**
 * @file NGPmcmc.h
 * @brief Non-Gaussian Process (NGP) MCMC sampler.
 *
 * Implements the Markov chain Monte Carlo (MCMC) algorithm for the locally
 * adaptive non-Gaussian process state-space model. Provides lifecycle functions,
 * matrix assembly routines, posterior parameter sampling, and file export helpers.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <lafp/gmlib.h>
#include <lafp/mdarray.h>
#include <lafp/SSsimulate2.h>

/**
 * @struct NGPmcmc
 * @brief Primary state container and workspace for NGP MCMC sampling.
 */
typedef struct {
    /* Inputs */
    gsl_matrix *y;        /**< Observation matrix (Nt x 1) */
    gsl_vector *tobs;     /**< Observation timestamps vector (Nt) */
    gsl_vector *delta;    /**< Time increment vector (Nt-1) */

    double sigEps;        /**< Observation noise standard deviation */
    gsl_vector *sigU;     /**< State process noise scale parameter */
    gsl_vector *sigA;     /**< Derivative process noise scale parameter */

    /* Sampler Configuration & Outputs */
    int Niter;            /**< Total MCMC sampling iterations */
    int burn;             /**< Burn-in iteration count */
    int Nt;               /**< Time series length */
    marray3d *th;         /**< Posterior state samples container (Niter x Nt x 3) */
    gsl_matrix *sig;      /**< Posterior scale samples container (Niter x 3) */

    /* Hyperparameter Priors */
    double a;             /**< Inverse-Gamma prior shape parameter */
    double b;             /**< Inverse-Gamma prior scale parameter */
    double sigMu;         /**< Prior standard deviation for mean state components */
    double sigAlph;       /**< Prior standard deviation for derivative state components */

    /* Temporary Workspaces & Sub-Samplers */
    gsl_vector *y_vec;

    /* State Simulation Workspace */
    marray3d *Z_s;
    marray3d *H_s;
    marray3d *T_s;
    marray3d *R_s;
    marray3d *Q_s;
    gsl_vector *a0_s;
    gsl_matrix *P0_s;
    SSsimulate2 *drawState;

    /* Observation Noise Workspace */
    gsl_vector *U;

    /* Acceptance Probability Workspace */
    double prob;
    gsl_vector *v_t;
    gsl_vector *vstar_t;
    gsl_vector *v_03;
    gsl_vector *v_02;
    gsl_vector *work3;
    gsl_vector *work2;
    gsl_vector *theta_t;
    gsl_vector *thetastar_t;
    gsl_vector *Gthetastar;
    gsl_matrix *HtildeT;

    /* Scale Sampler Workspace */
    gsl_vector *sigUstar;
    gsl_vector *sigAstar;
    gsl_matrix *theta;
    gsl_matrix *thetastar;
    gsl_vector *Ustar;
    gsl_vector *DUstar;
    gsl_vector *Astar;

    /* Matrix Assembly Dummies */
    marray3d *Zdummy;
    marray3d *Hdummy;
    marray3d *G;
    marray3d *H;
    marray3d *W;
    gsl_vector *a0dummy;
    gsl_matrix *P0dummy;

    marray3d *Ztdummy;
    marray3d *Htdummy;
    marray3d *Gtilde;
    marray3d *Htilde;
    marray3d *Wtilde;
    gsl_vector *a0tildedummy;
    gsl_matrix *P0tildedummy;

    marray3d *Zstardummy;
    marray3d *Hstardummy;
    marray3d *Gstardummy;
    marray3d *Tstardummy;
    marray3d *Wstar;
    gsl_vector *a0stardummy;
    gsl_matrix *P0stardummy;

    marray3d *Ztildestardummy;
    marray3d *Htildestardummy;
    marray3d *Gtildestardummy;
    marray3d *Ttildestardummy;
    marray3d *Wtildestar;
    gsl_vector *a0tildestardummy;
    gsl_matrix *P0tildestardummy;

    gsl_vector *a0_star;
    gsl_matrix *P0_star;
    marray3d *Z_star;
    marray3d *H_star;
    marray3d *T_star;
    marray3d *R_star;
    marray3d *Q_star;
    SSsimulate2 *drawThetaStar;

    gsl_rng *rand;        /**< GSL random number generator instance */

} NGPmcmc;

/***************************************
 * Constructor, destructor, main method *
 ***************************************/

/**
 * @brief Allocates a new NGPmcmc instance.
 * @return Pointer to allocated NGPmcmc object, or NULL on failure.
 */
NGPmcmc * NGPmcmc_New(void);

/**
 * @brief Initializes fields of an NGPmcmc instance to zero/NULL.
 * @param self Pointer to target NGPmcmc object.
 * @return 0 on success, non-zero on error.
 */
int NGPmcmc_init(NGPmcmc *self);

/**
 * @brief Frees all allocated memory within an NGPmcmc object and frees the struct itself.
 * @param s Pointer to target NGPmcmc object.
 * @return 0 on success, non-zero on error.
 */
int NGPmcmc_free(NGPmcmc *s);

/**
 * @brief Configures NGPmcmc sampler with observation data, iterations, and hyperparameters.
 * @param self Sampler instance pointer.
 * @param y Observation matrix (Nt x 1).
 * @param tobs Timestamps vector (Nt).
 * @param Niter Total sampling iterations.
 * @param sigU Process noise scale for states.
 * @param sigA Process noise scale for derivatives.
 * @param sigEps Initial observation noise standard deviation.
 * @param sigMu Prior scale for mean state components.
 * @param sigAlph Prior scale for derivative state components.
 * @param a Inverse-Gamma prior shape parameter.
 * @param b Inverse-Gamma prior scale parameter.
 * @param th Output 3D array for state samples (Niter x Nt x 3).
 * @param sig Output matrix for scale parameter samples (Niter x 3).
 * @return 0 on success, non-zero on error.
 */
int NGPmcmc_construct(NGPmcmc *self, gsl_matrix *y, gsl_vector *tobs, int Niter,
                      gsl_vector *sigU, gsl_vector *sigA, double sigEps,
                      double sigMu, double sigAlph, double a, double b,
                      marray3d *th, gsl_matrix *sig);

/**
 * @brief Runs the complete MCMC sampling loop over burn-in and saved iterations.
 * @param self Configured NGPmcmc sampler instance.
 * @return 0 on success, non-zero on error.
 */
int NGPmcmc_operations(NGPmcmc *self);

/*************************************
 * Sampling methods within operations *
 *************************************/

/** @brief Draws state vectors via state space smoothing simulation. */
int NGPmcmc_drawState(NGPmcmc *self);

/** @brief Draws observation noise variance sigEps from Inverse-Gamma conditional posterior. */
int NGPmcmc_drawSigEps(NGPmcmc *self);

/** @brief Computes Metropolis-Hastings acceptance probability for process noise candidate. */
int NGPmcmc_calc_acceptance_prob(NGPmcmc *self);

/** @brief Proposes and updates process noise scale parameters sigU and sigA. */
int NGPmcmc_drawSigs(NGPmcmc *self);

/****************************************
 * Convenience top-level execution entry *
 ****************************************/

/**
 * @brief Convenience one-shot function to allocate, run, export, and free NGP MCMC sampling.
 */
int NGPmcmc_NGPmcmc(gsl_matrix *y, gsl_vector *tobs, int Niter, gsl_vector *sigU,
                    gsl_vector *sigA, double sigEps, double sigMu, double sigAlph,
                    double a, double b, marray3d *th, gsl_matrix *sig);

/*********************************
 * Parameter modification methods *
 *********************************/

int NGPmcmc_updatePrior(NGPmcmc *self, char parameter, gsl_vector *param);

/*********************************
 * Output file export methods     *
 *********************************/

int NGPmcmc_writeTheta(NGPmcmc *self, char *name);
int NGPmcmc_writeSig(NGPmcmc *self, char *name);
int NGPmcmc_writeOutputs(NGPmcmc *self, char *baseName);

#endif
