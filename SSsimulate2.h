#ifndef _ABT_SSSIMULATE2_MVNORM_H
#define _ABT_SSSIMULATE2_MVNORM_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mdarray.h"
#include "fastStateSmoother2.h"
#include "KalmanFilter2.h"

typedef struct{
	Kalman2 *kalman1;
	FSS2 *fss1;
} KFS;

typedef struct{
	int Nt;
	int Np;
	int Nm;
	int Nr;
	int interleaved;

	//Vectors
	//Input
	gsl_vector *a_init;

	//Temporary
	gsl_vector *eps;
	gsl_vector *alphaplus0;
	gsl_vector *at_Vm;
	gsl_vector *Zat_Vp;
	gsl_vector *zerosNp;
	gsl_vector *zerosNr;
	gsl_vector *eta_Vr;
	gsl_vector *TTat_Vm;
	gsl_vector *RReta_Vm;

	//Matrices
	//Input
	gsl_matrix *y;
	gsl_matrix *P_init;
	gsl_matrix *alpha_draw;

	//Temporary
//	gsl_matrix *TT;
//	gsl_matrix *RR;
//	gsl_matrix *QQ;
//	gsl_matrix *HH;
//	gsl_matrix *ZZ;
	gsl_matrix *y_plus_Mtp;
	gsl_matrix *alpha_plus;
	gsl_matrix *alpha_hat;
	gsl_matrix *alpha_hat_plus;
	gsl_matrix *v;
	gsl_matrix *P_init_C;
	gsl_matrix *HH_C;
	gsl_matrix *QQ_C;

	//Marrays
	//Input
	marray3d *H;
	marray3d *T;
	marray3d *R;
	marray3d *Q;
	marray3d *Z;

	//Temporary
	marray3d *K1;
	marray3d *Finv1;
	marray3d *K2;
	marray3d *Finv2;

	//Random number generator
	gsl_rng *rand;

	//Kalman filters and fast state smoother
//	KFS *kfs1;
//	KFS *kfs2;
	Kalman2 *kalman1;
	Kalman2 *kalman2;
	FSS2 *fss1;
	FSS2 *fss2;

} SSsimulate2;

/*
*
*/

KFS * KFS_New();
int KFS_init(KFS *self);
int KFS_construct(KFS *self,gsl_matrix *y,gsl_vector *a_init,gsl_matrix *P_init,
                    marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q, 
                    gsl_matrix *v,marray3d *K,marray3d *Finv1,gsl_matrix *alpha_hat);
int KFS_free(KFS *self);
int KFS_operations(KFS *self);

/***************************************
* Constructor, destructor, main method *
***************************************/

SSsimulate2 * SSsimulate2_New();
int SSsimulate2_init(SSsimulate2 *self);
int SSsimulate2_construct(SSsimulate2 *self,gsl_matrix *y,gsl_vector *a_init,gsl_matrix *P_init,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,gsl_matrix *alpha_draw);
int SSsimulate2_free(SSsimulate2 *self);
int SSsimulate2_operations(SSsimulate2 *self);

/****************************************
* Class method for doing operation once *
****************************************/

int SSsimulate2_Simulate(gsl_matrix *y,gsl_vector *a_init,gsl_matrix *P_init,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,gsl_matrix *alpha_draw);

/*********************************
* Writing method outputs to file *
*********************************/

int SSsimulate2_writeAlphaDraw(SSsimulate2 *self,char *name);
int SSsimulate2_writeOutputs(SSsimulate2 *self,char *baseName);


#endif
