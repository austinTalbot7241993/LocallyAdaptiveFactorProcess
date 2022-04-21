#ifndef _ABT_KALMANFILTER2_H
#define _ABT_KALMANFILTER2_H

#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mdarray.h"

typedef struct{
	int Nt;
	int Np;
	int Nm;
	int Nr;

	//vectors
	gsl_vector *a_init;
	gsl_vector *yt_Vp;
	gsl_vector *at_Vm;
	gsl_vector *zat_Vp;
	gsl_vector *TTat_Vm;
	gsl_vector *KtVt_Vm;

	//matrices
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

	//arrays
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
* Constructor, destroctor, main method *
***************************************/

Kalman2 * Kalman2_New();
int Kalman2_init(Kalman2 *self);
int Kalman2_free(Kalman2 *s);
int Kalman2_construct(Kalman2 *self,gsl_matrix *y,gsl_vector *a_init,gsl_matrix *P_init,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,gsl_matrix *v,marray3d *K,marray3d *Finv);
int Kalman2_operations(Kalman2 *self);

/****************************************
* Class method for doing operation once *
****************************************/

int Kalman2_Filter(gsl_matrix *y,gsl_vector *a_init,gsl_matrix *P_init,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,gsl_matrix *v,marray3d *K,marray3d *Finv);

/*********************************
* Writing method outputs to file *
*********************************/

int Kalman2_writeV(Kalman2 *self,char *name);
int Kalman2_writeK(Kalman2 *self,char *name);
int Kalman2_writeFinv(Kalman2 *self,char *name);
int Kalman2_writeOutputs(Kalman2 *self,char *baseName);


#endif
