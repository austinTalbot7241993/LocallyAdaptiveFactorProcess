#ifndef _ABT_KALMANFILTER2_H
#define _ABT_KALMANFILTER2_H

/**
 * @file KalmanFilter2.h
 * @brief Forward Kalman Filtering module for state-space time series models.
 *
 * Implements forward Kalman filtering recursions as described in Durbin & Koopman (2012).
 */

#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lafp/mdarray.h>

/**
 * @struct Kalman2
 * @brief Forward Kalman Filter state container and scratch workspace.
 */
typedef struct {
    int Nt;               /**< Time series length */
    int Np;               /**< Observation dimension */
    int Nm;               /**< State vector dimension */
    int Nr;               /**< State noise vector dimension */

    /* Vectors */
    gsl_vector *a_init;
    gsl_vector *yt_Vp;
    gsl_vector *at_Vm;
    gsl_vector *zat_Vp;
    gsl_vector *TTat_Vm;
    gsl_vector *KtVt_Vm;

    /* Matrices */
    gsl_matrix *a_Mtm;
    gsl_matrix *P_init;
    gsl_matrix *v;
    gsl_matrix *y;
    gsl_matrix *ZPt_Mpm;
    gsl_matrix *Zt_Mmp;
    gsl_matrix *ZPZ_Mpp;
    gsl_matrix *Finvt_Mpp;
    gsl_matrix *TTPtZt_Mmp;
    gsl_matrix *TTPt_Mmm;
    gsl_matrix *Kt_Mmp;
    gsl_matrix *L_Mmm;
    gsl_matrix *KtZ_Mmm;
    gsl_matrix *Lt_Mmm;
    gsl_matrix *TTPtLt_Mmm;
    gsl_matrix *RRQQ;
    gsl_matrix *RQR;
    gsl_matrix *RRT;

    /* 3D Arrays */
    marray3d *P_Atmm;
    marray3d *T;
    marray3d *R;
    marray3d *Q;
    marray3d *K;
    marray3d *Finv;
    marray3d *H;
    marray3d *Z;

} Kalman2;

/***************************************
 * Constructor, destructor, main method *
 ***************************************/

/**
 * Instantiate an uninitialized Kalman2 filter object.
 *
 * Returns
 * -------
 * Kalman2 *
 *     Pointer to allocated Kalman2 object, or NULL on allocation failure.
 */
Kalman2 * Kalman2_New(void);

/**
 * Initialize fields of a Kalman2 struct to zero/NULL.
 *
 * Parameters
 * ----------
 * self : Kalman2 *
 *     Target Kalman2 instance.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int Kalman2_init(Kalman2 *self);

/**
 * Free all internal matrix workspaces and deallocate the Kalman2 struct.
 *
 * Parameters
 * ----------
 * s : Kalman2 *
 *     Target Kalman2 instance pointer to release.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int Kalman2_free(Kalman2 *s);

/**
 * Configure and pre-allocate workspaces for forward Kalman filtering.
 *
 * Parameters
 * ----------
 * self : Kalman2 *
 *     Target Kalman2 filter instance.
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
 * v : gsl_matrix *
 *     Output innovation residuals matrix of shape (Nt, Np).
 * K : marray3d *
 *     Output Kalman gain matrices array of shape (Nt, Nm, Np).
 * Finv : marray3d *
 *     Output prediction covariance inverse matrices array of shape (Nt, Np, Np).
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int Kalman2_construct(Kalman2 *self, gsl_matrix *y, gsl_vector *a_init,
                      gsl_matrix *P_init, marray3d *Z, marray3d *H,
                      marray3d *T, marray3d *R, marray3d *Q,
                      gsl_matrix *v, marray3d *K, marray3d *Finv);

/**
 * Execute forward Kalman filtering recursions over time steps t = 0 ... Nt-1.
 *
 * Parameters
 * ----------
 * self : Kalman2 *
 *     Configured Kalman2 filter instance.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 *
 * Notes
 * -----
 * For each time step t = 0 ... Nt-1, computes:
 * 1. Residuals: v_t = y_t - Z_t * a_t
 * 2. Forecast Covariance: F_t = Z_t * P_t * Z_t^T + H_t
 * 3. Kalman Gain: K_t = T_t * P_t * Z_t^T * F_t^{-1}
 * 4. State Update: a_{t+1} = T_t * a_t + K_t * v_t
 * 5. Covariance Update: P_{t+1} = T_t * P_t * L_t^T + R_t * Q_t * R_t^T
 */
int Kalman2_operations(Kalman2 *self);

/****************************************
 * Convenience top-level execution entry *
 ****************************************/

/**
 * One-shot function to allocate, run forward Kalman filter, and free resources.
 *
 * Parameters
 * ----------
 * y : gsl_matrix *
 *     Observation matrix of shape (Nt, Np).
 * a_init : gsl_vector *
 *     Initial state mean of shape (Nm,).
 * P_init : gsl_matrix *
 *     Initial state covariance of shape (Nm, Nm).
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
 * v : gsl_matrix *
 *     Output innovation matrix of shape (Nt, Np).
 * K : marray3d *
 *     Output Kalman gain array.
 * Finv : marray3d *
 *     Output forecast covariance inverse array.
 *
 * Returns
 * -------
 * int
 *     0 if successful, non-zero error code on failure.
 */
int Kalman2_Filter(gsl_matrix *y, gsl_vector *a_init, gsl_matrix *P_init,
                   marray3d *Z, marray3d *H, marray3d *T, marray3d *R,
                   marray3d *Q, gsl_matrix *v, marray3d *K, marray3d *Finv);

/*********************************
 * Writing method outputs to file *
 *********************************/

int Kalman2_writeV(Kalman2 *self, char *name);
int Kalman2_writeK(Kalman2 *self, char *name);
int Kalman2_writeFinv(Kalman2 *self, char *name);
int Kalman2_writeOutputs(Kalman2 *self, char *baseName);

#endif
