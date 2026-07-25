#ifndef _ABT_NGPMCMC_H
#define _ABT_NGPMCMC_H

/**
 * @file NGPmcmc.h
 * @brief Non-Gaussian Process (NGP) MCMC sampler module.
 *
 * Implements the Markov chain Monte Carlo (MCMC) algorithm for the locally
 * adaptive non-Gaussian process state-space model described in Durante & Dunson (2014).
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
 * Instantiate an uninitialized NGPmcmc sampler.
 *
 * Returns
 * -------
 * NGPmcmc *
 *     Pointer to newly allocated NGPmcmc instance, or NULL on allocation failure.
 */
NGPmcmc * NGPmcmc_New(void);

/**
 * Initialize fields of an NGPmcmc sampler struct to zero/NULL.
 *
 * Parameters
 * ----------
 * self : NGPmcmc *
 *     Target NGPmcmc instance pointer.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int NGPmcmc_init(NGPmcmc *self);

/**
 * Free all internal matrix workspaces and deallocate the NGPmcmc sampler struct.
 *
 * Parameters
 * ----------
 * s : NGPmcmc *
 *     Target NGPmcmc instance pointer to release.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int NGPmcmc_free(NGPmcmc *s);

/**
 * Configure and pre-allocate workspaces for NGP MCMC model fitting.
 *
 * Parameters
 * ----------
 * self : NGPmcmc *
 *     Target NGPmcmc sampler instance.
 * y : gsl_matrix *
 *     Observation matrix of shape (Nt, 1) or (Nt, Np).
 * tobs : gsl_vector *
 *     Timestamps vector of shape (Nt,).
 * Niter : int
 *     Total number of MCMC sampling iterations to execute.
 * sigU : gsl_vector *
 *     Initial process scale parameter vector for state components.
 * sigA : gsl_vector *
 *     Initial process scale parameter vector for derivative components.
 * sigEps : double
 *     Initial observation noise standard deviation.
 * sigMu : double
 *     Prior scale standard deviation for mean state components.
 * sigAlph : double
 *     Prior scale standard deviation for derivative state components.
 * a : double
 *     Inverse-Gamma prior shape parameter for noise variance.
 * b : double
 *     Inverse-Gamma prior scale parameter for noise variance.
 * th : marray3d *
 *     Output 3D array of shape (Niter, Nt, 3) storing posterior state draws.
 * sig : gsl_matrix *
 *     Output matrix of shape (Niter, 3) storing posterior scale draws.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 *
 * Notes
 * -----
 * Memory allocations for internal matrix workspaces occur exclusively inside
 * NGPmcmc_construct. No dynamic memory allocation or free occurs during the
 * main MCMC loop (NGPmcmc_operations).
 */
int NGPmcmc_construct(NGPmcmc *self, gsl_matrix *y, gsl_vector *tobs, int Niter,
                      gsl_vector *sigU, gsl_vector *sigA, double sigEps,
                      double sigMu, double sigAlph, double a, double b,
                      marray3d *th, gsl_matrix *sig);

/**
 * Run the complete MCMC sampling loop over burn-in and saved iterations.
 *
 * Parameters
 * ----------
 * self : NGPmcmc *
 *     Configured NGPmcmc sampler instance.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 *
 * Notes
 * -----
 * Executes Niter sampling loops. Each iteration performs sequential draws:
 * 1. NGPmcmc_drawState() - Durbin-Koopman simulation smoother.
 * 2. NGPmcmc_drawSigEps() - Inverse-Gamma draw for observation variance.
 * 3. NGPmcmc_drawSigs() - Metropolis-Hastings update for process noise scales.
 */
int NGPmcmc_operations(NGPmcmc *self);

/*************************************
 * Sampling methods within operations *
 *************************************/

/**
 * Draw state trajectory theta from conditional posterior via simulation smoother.
 *
 * Parameters
 * ----------
 * self : NGPmcmc *
 *     Configured NGPmcmc sampler instance.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int NGPmcmc_drawState(NGPmcmc *self);

/**
 * Draw observation noise standard deviation sigEps from Inverse-Gamma posterior.
 *
 * Parameters
 * ----------
 * self : NGPmcmc *
 *     Configured NGPmcmc sampler instance.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int NGPmcmc_drawSigEps(NGPmcmc *self);

/**
 * Compute Metropolis-Hastings acceptance probability for proposed process noise candidate.
 *
 * Parameters
 * ----------
 * self : NGPmcmc *
 *     Configured NGPmcmc sampler instance.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int NGPmcmc_calc_acceptance_prob(NGPmcmc *self);

/**
 * Propose and update process noise scale parameters sigU and sigA via Metropolis-Hastings.
 *
 * Parameters
 * ----------
 * self : NGPmcmc *
 *     Configured NGPmcmc sampler instance.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int NGPmcmc_drawSigs(NGPmcmc *self);

/****************************************
 * Convenience top-level execution entry *
 ****************************************/

/**
 * One-shot function to allocate, configure, run MCMC, export results, and free resources.
 *
 * Parameters
 * ----------
 * y : gsl_matrix *
 *     Observation matrix of shape (Nt, 1).
 * tobs : gsl_vector *
 *     Timestamps vector of shape (Nt,).
 * Niter : int
 *     Total number of MCMC sampling iterations.
 * sigU : gsl_vector *
 *     Initial process scale for states.
 * sigA : gsl_vector *
 *     Initial process scale for derivatives.
 * sigEps : double
 *     Initial observation noise standard deviation.
 * sigMu : double
 *     Prior scale standard deviation for mean state components.
 * sigAlph : double
 *     Prior scale standard deviation for derivative state components.
 * a : double
 *     Inverse-Gamma prior shape parameter.
 * b : double
 *     Inverse-Gamma prior scale parameter.
 * th : marray3d *
 *     Output 3D array of shape (Niter, Nt, 3) storing posterior state draws.
 * sig : gsl_matrix *
 *     Output matrix of shape (Niter, 3) storing posterior scale draws.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
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

/**
 * Export MCMC posterior draws to text files on disk.
 *
 * Parameters
 * ----------
 * self : NGPmcmc *
 *     Sampler instance containing posterior draws.
 * baseName : char *
 *     Base filename prefix (e.g. "my_results" produces "my_results_Theta.txt" and "my_results_Sig.txt").
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int NGPmcmc_writeOutputs(NGPmcmc *self, char *baseName);

#endif
