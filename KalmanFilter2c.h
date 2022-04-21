#ifndef _ABT_KALMANFILTER2C_H
#define _ABT_KALMANFILTER2C_H

#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mdarray_complex.h"

typedef struct{
	int Nt;
	int Np;
	int Nm;
	int Nr;

	//vectors
	gsl_vector_complex *a_init;
	gsl_vector_complex *yt_Vp;
	gsl_vector_complex *at_Vm;
	gsl_vector_complex *zat_Vp;
	gsl_vector_complex *TTat_Vm;
	gsl_vector_complex *KtVt_Vm;

	//matrices
	gsl_matrix_complex *a_Mtm;
	gsl_matrix_complex *P_init;
	gsl_matrix_complex *v;
	gsl_matrix_complex *y;
	gsl_matrix_complex *ZPt_Mpm;
	gsl_matrix_complex *Zt_Mmp;
	gsl_matrix_complex *ZPZ_Mpp;
	gsl_matrix_complex *Finvt_Mpp;
	gsl_matrix_complex *TTPtZt_Mmp;
	gsl_matrix_complex *TTPt_Mmm;
	gsl_matrix_complex *Kt_Mmp;
	gsl_matrix_complex *L_Mmm;
	gsl_matrix_complex *KtZ_Mmm;
	gsl_matrix_complex *Lt_Mmm;
	gsl_matrix_complex *TTPtLt_Mmm;
	gsl_matrix_complex *RRQQ;
	gsl_matrix_complex *RQR;
	gsl_matrix_complex *RRT;

	//arrays
	marray3d_complex *P_Atmm;
	marray3d_complex *T;
	marray3d_complex *R;
	marray3d_complex *Q;
	marray3d_complex *K;
	marray3d_complex *Finv;
	marray3d_complex *H;
	marray3d_complex *Z;

} Kalman2c;

/***************************************
* Constructor, destroctor, main method *
***************************************/

Kalman2c * Kalman2c_New();
int Kalman2c_init(Kalman2c *self);
int Kalman2c_free(Kalman2c *s);
int Kalman2c_construct(Kalman2c *self,gsl_matrix_complex *y,gsl_vector_complex *a_init,gsl_matrix_complex *P_init,marray3d_complex *Z,marray3d_complex *H,marray3d_complex *T,marray3d_complex *R,marray3d_complex *Q,gsl_matrix_complex *v,marray3d_complex *K,marray3d_complex *Finv);
int Kalman2c_operations(Kalman2c *self);

/****************************************
* Class method for doing operation once *
****************************************/

int Kalman2c_Filter(gsl_matrix_complex *y,gsl_vector_complex *a_init,gsl_matrix_complex *P_init,marray3d_complex *Z,marray3d_complex *H,marray3d_complex *T,marray3d_complex *R,marray3d_complex *Q,gsl_matrix_complex *v,marray3d_complex *K,marray3d_complex *Finv);

/*********************************
* Writing method outputs to file *
*********************************/

int Kalman2c_writeV(Kalman2c *self,char *name);
int Kalman2c_writeK(Kalman2c *self,char *name);
int Kalman2c_writeFinv(Kalman2c *self,char *name);
int Kalman2c_writeOutputs(Kalman2c *self,char *baseName);


#endif
