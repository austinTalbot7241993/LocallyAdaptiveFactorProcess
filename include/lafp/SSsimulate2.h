#ifndef _ABT_SSSIMULATE2_MVNORM_H
#define _ABT_SSSIMULATE2_MVNORM_H

/**
 * @file SSsimulate2.h
 * @brief State-Space Simulation Smoother module (Durbin & Koopman simulation smoother).
 *
 * Implements the simulation smoother algorithm for state-space models using paired
 * Kalman filtering and Fast State Smoothing (FSS).
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lafp/mdarray.h>
#include <lafp/fastStateSmoother2.h>
#include <lafp/KalmanFilter2.h>

typedef struct {
    Kalman2 *kalman1;
    FSS2 *fss1;
} KFS;

/**
 * @struct SSsimulate2
 * @brief Primary workspace container for state-space simulation smoothing.
 */
typedef struct {
    int Nt;               /**< Time series length */
    int Np;               /**< Observation dimension */
    int Nm;               /**< State vector dimension */
    int Nr;               /**< State noise vector dimension */
    int interleaved;      /**< Interleaved execution flag */

    /* Input Vectors */
    gsl_vector *a_init;

    /* Temporary Vectors */
    gsl_vector *eps;
    gsl_vector *alphaplus0;
    gsl_vector *at_Vm;
    gsl_vector *Zat_Vp;
    gsl_vector *zerosNp;
    gsl_vector *zerosNr;
    gsl_vector *eta_Vr;
    gsl_vector *TTat_Vm;
    gsl_vector *RReta_Vm;

    /* Input Matrices */
    gsl_matrix *y;
    gsl_matrix *P_init;
    gsl_matrix *alpha_draw;

    /* Temporary Matrices */
    gsl_matrix *y_plus_Mtp;
    gsl_matrix *alpha_plus;
    gsl_matrix *alpha_hat;
    gsl_matrix *alpha_hat_plus;
    gsl_matrix *v;
    gsl_matrix *P_init_C;
    gsl_matrix *HH_C;
    gsl_matrix *QQ_C;

    /* Input 3D Arrays */
    marray3d *H;
    marray3d *T;
    marray3d *R;
    marray3d *Q;
    marray3d *Z;

    /* Temporary 3D Arrays */
    marray3d *K1;
    marray3d *Finv1;
    marray3d *K2;
    marray3d *Finv2;

    gsl_rng *rand;        /**< GSL random number generator instance */

    /* Paired Kalman Filters & Smoothers */
    Kalman2 *kalman1;
    Kalman2 *kalman2;
    FSS2 *fss1;
    FSS2 *fss2;

} SSsimulate2;

/* KFS Helper Methods */
KFS * KFS_New(void);
int KFS_init(KFS *self);
int KFS_construct(KFS *self, gsl_matrix *y, gsl_vector *a_init, gsl_matrix *P_init,
                  marray3d *Z, marray3d *H, marray3d *T, marray3d *R, marray3d *Q,
                  gsl_matrix *v, marray3d *K, marray3d *Finv1, gsl_matrix *alpha_hat);
int KFS_free(KFS *self);
int KFS_operations(KFS *self);

/***************************************
 * Constructor, destructor, main method *
 ***************************************/

/**
 * Instantiate an uninitialized SSsimulate2 simulation smoother instance.
 *
 * Returns
 * -------
 * SSsimulate2 *
 *     Pointer to newly allocated SSsimulate2 instance, or NULL on allocation failure.
 */
SSsimulate2 * SSsimulate2_New(void);

/**
 * Initialize fields of an SSsimulate2 instance to zero/NULL.
 *
 * Parameters
 * ----------
 * self : SSsimulate2 *
 *     Target SSsimulate2 instance.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int SSsimulate2_init(SSsimulate2 *self);

/**
 * Configure and pre-allocate workspace for state-space simulation smoothing.
 *
 * Parameters
 * ----------
 * self : SSsimulate2 *
 *     Target simulation smoother instance.
 * y : gsl_matrix *
 *     Observation matrix of shape (Nt, Np).
 * a_init : gsl_vector *
 *     Initial state prior mean vector of shape (Nm,).
 * P_init : gsl_matrix *
 *     Initial state prior covariance matrix of shape (Nm, Nm).
 * Z : marray3d *
 *     Observation design matrix array of shape (Nt, Np, Nm) or (1, Np, Nm).
 * H : marray3d *
 *     Observation noise covariance array of shape (Nt, Np, Np) or (1, Np, Np).
 * T : marray3d *
 *     State transition matrix array of shape (Nt, Nm, Nm) or (1, Nm, Nm).
 * R : marray3d *
 *     State disturbance selection matrix array of shape (Nt, Nm, Nr) or (1, Nm, Nr).
 * Q : marray3d *
 *     State disturbance covariance matrix array of shape (Nt, Nr, Nr) or (1, Nr, Nr).
 * alpha_draw : gsl_matrix *
 *     Output matrix of shape (Nt, Nm) where the drawn state trajectory will be written.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int SSsimulate2_construct(SSsimulate2 *self, gsl_matrix *y, gsl_vector *a_init,
                          gsl_matrix *P_init, marray3d *Z, marray3d *H,
                          marray3d *T, marray3d *R, marray3d *Q,
                          gsl_matrix *alpha_draw);

/**
 * Free all internal matrix workspaces and deallocate the SSsimulate2 struct.
 *
 * Parameters
 * ----------
 * self : SSsimulate2 *
 *     Target simulation smoother instance pointer to release.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int SSsimulate2_free(SSsimulate2 *self);

/**
 * Run Durbin-Koopman simulation smoother pass to sample state trajectory.
 *
 * Parameters
 * ----------
 * self : SSsimulate2 *
 *     Configured SSsimulate2 instance pointer.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 *
 * Notes
 * -----
 * Executes two paired passes of Kalman filtering and Fast State Smoothing (FSS):
 * 1. Simulates pseudo-random state trajectories alpha+ and observations y+.
 * 2. Filters y+ and original observations y using dual Kalman filters.
 * 3. Smooths states via dual Fast State Smoothers to construct exact draw alpha_draw.
 */
int SSsimulate2_operations(SSsimulate2 *self);

/****************************************
 * Convenience top-level simulation entry *
 ****************************************/

/**
 * One-shot function to allocate, run simulation smoother, and free resources.
 *
 * Parameters
 * ----------
 * y : gsl_matrix *
 *     Observation matrix of shape (Nt, Np).
 * a_init : gsl_vector *
 *     Initial state prior mean of shape (Nm,).
 * P_init : gsl_matrix *
 *     Initial state prior covariance of shape (Nm, Nm).
 * Z : marray3d *
 *     Observation design matrix array.
 * H : marray3d *
 *     Observation noise covariance array.
 * T : marray3d *
 *     State transition matrix array.
 * R : marray3d *
 *     State disturbance selection matrix array.
 * Q : marray3d *
 *     State disturbance covariance array.
 * alpha_draw : gsl_matrix *
 *     Output state trajectory matrix of shape (Nt, Nm).
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int SSsimulate2_Simulate(gsl_matrix *y, gsl_vector *a_init, gsl_matrix *P_init,
                         marray3d *Z, marray3d *H, marray3d *T, marray3d *R,
                         marray3d *Q, gsl_matrix *alpha_draw);

/*********************************
 * File export methods            *
 *********************************/

int SSsimulate2_writeAlphaDraw(SSsimulate2 *self, char *name);
int SSsimulate2_writeOutputs(SSsimulate2 *self, char *baseName);

#endif
