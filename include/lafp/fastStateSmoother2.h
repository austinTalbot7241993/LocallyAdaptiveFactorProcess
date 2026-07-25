#ifndef _ABT_FASTSTATESMOOTHER2_H
#define _ABT_FASTSTATESMOOTHER2_H

/**
 * @file fastStateSmoother2.h
 * @brief Fast State Smoothing (FSS) backward recursions for state-space models.
 *
 * Implements the linear-time Fast State Smoother (Durbin & Koopman 2012) for computing
 * smoothed state vector estimates and disturbances.
 */

#include <gsl/gsl_matrix.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <lafp/mdarray.h>

/**
 * @struct FSS2
 * @brief Fast State Smoother workspace container.
 */
typedef struct {
    int Nt;               /**< Time series length */
    int Np;               /**< Observation dimension */
    int Nm;               /**< State vector dimension */
    int Nr;               /**< State noise vector dimension */

    /* Input Vectors */
    gsl_vector *a_init;

    /* Temporary Vectors */
    gsl_vector *vt_Vp;
    gsl_vector *Fvt_Vp;
    gsl_vector *Kttrt_Vp;
    gsl_vector *rt_Vm;
    gsl_vector *r_init;
    gsl_vector *TTrt_Vm;
    gsl_vector *ZZFvt_Vm;
    gsl_vector *Pr_Vm;
    gsl_vector *TTat_Vm;
    gsl_vector *RQRrt_Vm;
    gsl_vector *at_Vm;

    /* Input Matrices */
    gsl_matrix *v;
    gsl_matrix *P_init;
    gsl_matrix *alpha;

    /* Temporary Matrices */
    gsl_matrix *r;
    gsl_matrix *ZZT_Mmp;
    gsl_matrix *RRQQ;
    gsl_matrix *RRt;
    gsl_matrix *RQR;
    gsl_matrix *Ktt_Mpm;
    gsl_matrix *TTt_Mmm;

    /* Input 3D Arrays */
    marray3d *Z;
    marray3d *T;
    marray3d *R;
    marray3d *Q;
    marray3d *K;
    marray3d *Finv;
} FSS2;

/***************************************
 * Constructor, destructor, main method *
 ***************************************/

/**
 * Instantiate an uninitialized FSS2 smoother object.
 *
 * Returns
 * -------
 * FSS2 *
 *     Pointer to allocated FSS2 object, or NULL on allocation failure.
 */
FSS2 * FSS2_New(void);

/**
 * Initialize fields of an FSS2 struct to zero/NULL.
 *
 * Parameters
 * ----------
 * self : FSS2 *
 *     Target FSS2 instance.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int FSS2_init(FSS2 *self);

/**
 * Free all internal matrix workspaces and deallocate the FSS2 struct.
 *
 * Parameters
 * ----------
 * s : FSS2 *
 *     Target FSS2 instance pointer to release.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int FSS2_free(FSS2 *s);

/**
 * Configure and pre-allocate workspaces for Fast State Smoothing.
 *
 * Parameters
 * ----------
 * self : FSS2 *
 *     Target FSS2 smoother instance.
 * v : gsl_matrix *
 *     Innovation residuals matrix from Kalman filter of shape (Nt, Np).
 * K : marray3d *
 *     Kalman gain array from Kalman filter of shape (Nt, Nm, Np).
 * Finv : marray3d *
 *     Forecast covariance inverse array from Kalman filter of shape (Nt, Np, Np).
 * a_init : gsl_vector *
 *     Initial state mean vector of shape (Nm,).
 * P_init : gsl_matrix *
 *     Initial state covariance matrix of shape (Nm, Nm).
 * Z : marray3d *
 *     Observation design matrix array.
 * T : marray3d *
 *     State transition matrix array.
 * R : marray3d *
 *     State disturbance selection matrix array.
 * Q : marray3d *
 *     State disturbance covariance array.
 * alpha : gsl_matrix *
 *     Output matrix of shape (Nt, Nm) where smoothed states are stored.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int FSS2_construct(FSS2 *self, gsl_matrix *v, marray3d *K, marray3d *Finv,
                   gsl_vector *a_init, gsl_matrix *P_init, marray3d *Z,
                   marray3d *T, marray3d *R, marray3d *Q, gsl_matrix *alpha);

/**
 * Execute backward Fast State Smoothing recursions over t = Nt-1 ... 0.
 *
 * Parameters
 * ----------
 * self : FSS2 *
 *     Configured FSS2 smoother instance.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 *
 * Notes
 * -----
 * Backward recursions compute smoothing vector r_t and smoothed states alpha_t:
 * 1. u_t = F_t^{-1} * v_t - K_t^T * r_t
 * 2. r_{t-1} = Z_t^T * u_t + T_t^T * r_t
 * 3. \hat{\alpha}_1 = a_1 + P_1 * r_0
 * 4. \hat{\alpha}_{t+1} = T_t * \hat{\alpha}_t + R_t * Q_t * R_t^T * r_t
 */
int FSS2_operations(FSS2 *self);

/****************************************
 * Convenience top-level execution entry *
 ****************************************/

/**
 * One-shot function to allocate, run Fast State Smoother, and free resources.
 *
 * Parameters
 * ----------
 * v : gsl_matrix *
 *     Innovation residuals matrix of shape (Nt, Np).
 * K : marray3d *
 *     Kalman gain array.
 * Finv : marray3d *
 *     Forecast covariance inverse array.
 * a_init : gsl_vector *
 *     Initial state mean.
 * P_init : gsl_matrix *
 *     Initial state covariance.
 * Z : marray3d *
 *     Observation design matrix array.
 * T : marray3d *
 *     State transition matrix array.
 * R : marray3d *
 *     State disturbance selection matrix array.
 * Q : marray3d *
 *     State disturbance covariance array.
 * alpha : gsl_matrix *
 *     Output smoothed state matrix of shape (Nt, Nm).
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int FSS2_Smoother(gsl_matrix *v, marray3d *K, marray3d *Finv, gsl_vector *a_init,
                  gsl_matrix *P_init, marray3d *Z, marray3d *T, marray3d *R,
                  marray3d *Q, gsl_matrix *alpha);

/*********************************
 * Writing method outputs to file *
 *********************************/

int FSS2_writeAlpha(FSS2 *self, char *name);
int FSS2_writeOutputs(FSS2 *self, char *baseName);

#endif
