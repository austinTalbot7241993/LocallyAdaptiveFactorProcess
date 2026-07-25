#ifndef _ABT_SSSIMULATE2_MVNORM_H
#define _ABT_SSSIMULATE2_MVNORM_H

/**
 * @file SSsimulate2.h
 * @brief State-Space Simulation Smoother algorithm (Durbin & Koopman simulation smoother).
 *
 * Provides routines for sampling disturbance vectors and state vectors in state-space
 * models using paired Kalman filtering and Fast State Smoothing (FSS).
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lafp/mdarray.h>
#include <lafp/fastStateSmoother2.h>
#include <lafp/KalmanFilter2.h>

/**
 * @struct KFS
 * @brief Pair container for Kalman Filter and Fast State Smoother handles.
 */
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
 * @brief Allocates a new SSsimulate2 object.
 * @return Pointer to allocated instance or NULL.
 */
SSsimulate2 * SSsimulate2_New(void);

/**
 * @brief Initializes fields of an SSsimulate2 object to zero/NULL.
 */
int SSsimulate2_init(SSsimulate2 *self);

/**
 * @brief Configures SSsimulate2 instance with state-space model matrices.
 */
int SSsimulate2_construct(SSsimulate2 *self, gsl_matrix *y, gsl_vector *a_init,
                          gsl_matrix *P_init, marray3d *Z, marray3d *H,
                          marray3d *T, marray3d *R, marray3d *Q,
                          gsl_matrix *alpha_draw);

/**
 * @brief Frees all allocated memory in an SSsimulate2 instance and frees the struct.
 */
int SSsimulate2_free(SSsimulate2 *self);

/**
 * @brief Runs simulation smoother pass to sample state trajectory alpha_draw.
 */
int SSsimulate2_operations(SSsimulate2 *self);

/****************************************
 * Convenience top-level simulation entry *
 ****************************************/

/**
 * @brief One-shot function to allocate, run simulation smoother, and free resources.
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
