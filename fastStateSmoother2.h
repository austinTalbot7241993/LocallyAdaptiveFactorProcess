#ifndef _ABT_FASTSTATESMOOTHER2_H
#define _ABT_FASTSTATESMOOTHER2_H

#include <gsl/gsl_matrix.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mdarray.h"

typedef struct{
	int Nt;
	int Np;
	int Nm;
	int Nr;

	//Vectors
	//Input
	gsl_vector *a_init;

	//Temporary
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
	
	//Matrices
	//Input
	gsl_matrix *v;
	gsl_matrix *P_init;
	gsl_matrix *alpha;

	//Temporary
	gsl_matrix *r;
	gsl_matrix *ZZT_Mmp;
	gsl_matrix *RRQQ;
	gsl_matrix *RRt;
	gsl_matrix *RQR;
	gsl_matrix *Ktt_Mpm;
	gsl_matrix *TTt_Mmm;

	//Arrays
	//Input
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

FSS2 * FSS2_New();
int FSS2_init(FSS2 *self);
int FSS2_free(FSS2 *s);
int FSS2_construct(FSS2 *self,gsl_matrix *v,marray3d *K,marray3d *Finv,gsl_vector *a_init,gsl_matrix *P_init,marray3d *Z,marray3d *T,marray3d *R,marray3d *Q,gsl_matrix *alpha);
int FSS2_operations(FSS2 *self);

/****************************************
* Class method for doing operation once *
****************************************/

int FSS2_Smoother(gsl_matrix *v,marray3d *K,marray3d *Finv,gsl_vector *a_init,gsl_matrix *P_init,marray3d *Z,marray3d *T,marray3d *R,marray3d *Q,gsl_matrix *alpha);

/*********************************
* Writing method outputs to file *
*********************************/

int FSS2_writeAlpha(FSS2 *self,char *name);
int FSS2_writeOutputs(FSS2 *self,char *baseName);

#endif
