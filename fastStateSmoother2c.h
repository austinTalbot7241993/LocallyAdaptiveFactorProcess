#ifndef _ABT_FASTSTATESMOOTHER2C_H
#define _ABT_FASTSTATESMOOTHER2C_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_complex.h>
#include <gsl/gsl_complex.h>
#include "mdarray_complex.h"

typedef struct{
	int Nt;
	int Np;
	int Nm;
	int Nr;

	//Vectors
	//Input
	gsl_vector_complex *a_init;

	//Temporary
	gsl_vector_complex *vt_Vp;
	gsl_vector_complex *Fvt_Vp;
	gsl_vector_complex *Kttrt_Vp;
	gsl_vector_complex *rt_Vm;
	gsl_vector_complex *r_init;
	gsl_vector_complex *TTrt_Vm;
	gsl_vector_complex *ZZFvt_Vm;
	gsl_vector_complex *Pr_Vm;
	gsl_vector_complex *TTat_Vm;
	gsl_vector_complex *RQRrt_Vm;
	gsl_vector_complex *at_Vm;
	
	//Matrices
	//Input
	gsl_matrix_complex *v;
	gsl_matrix_complex *P_init;
	gsl_matrix_complex *alpha;

	//Temporary
	gsl_matrix_complex *r;
	gsl_matrix_complex *ZZT_Mmp;
	gsl_matrix_complex *RRQQ;
	gsl_matrix_complex *RRt;
	gsl_matrix_complex *RQR;
	gsl_matrix_complex *Ktt_Mpm;
	gsl_matrix_complex *TTt_Mmm;

	//Arrays
	//Input
	marray3d_complex *Z;
	marray3d_complex *T;
	marray3d_complex *R;
	marray3d_complex *Q;
	marray3d_complex *K;
	marray3d_complex *Finv;
} FSS2c;

/***************************************
* Constructor, destructor, main method *
***************************************/

FSS2c * FSS2c_New();
int FSS2c_init(FSS2c *self);
int FSS2c_free(FSS2c *s);
int FSS2c_construct(FSS2c *self,gsl_matrix_complex *v,marray3d_complex *K,marray3d_complex *Finv,gsl_vector_complex *a_init,gsl_matrix_complex *P_init,marray3d_complex *Z,marray3d_complex *T,marray3d_complex *R,marray3d_complex *Q,gsl_matrix_complex *alpha);
int FSS2c_operations(FSS2c *self);

/****************************************
* Class method for doing operation once *
****************************************/

int FSS2c_Smoother(gsl_matrix_complex *v,marray3d_complex *K,marray3d_complex *Finv,gsl_vector_complex *a_init,gsl_matrix_complex *P_init,marray3d_complex *Z,marray3d_complex *T,marray3d_complex *R,marray3d_complex *Q,gsl_matrix_complex *alpha);

/*********************************
* Writing method outputs to file *
*********************************/

int FSS2c_writeAlpha(FSS2c *self,char *name);
int FSS2c_writeOutputs(FSS2c *self,char *baseName);

#endif
