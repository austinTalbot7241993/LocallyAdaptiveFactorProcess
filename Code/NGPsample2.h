#ifndef _ABT_NGPSAMPLE2_H
#define _ABT_NGPSAMPLE2_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include "mdarray.h"

typedef struct{
	int Np;
	int Ns;
	int Nt;
	int Nm;
	int Nr;

	gsl_vector *a0;

	gsl_matrix *y;
	gsl_matrix *P0;
	gsl_matrix *alphaDraw;

	marray3d *Z;
	marray3d *T;
	marray3d *R;
	marray3d *H;
	marray3d *Q;
	marray3d *alphaSamples;

	SSsimulate2 *SSsample;

} _NGPsample;

typedef struct{
	int Nr;
	int Np;
	int Ns;
	int Nt;
	int approx;
	int Nsamples;

	gsl_vector *a0;

	gsl_matrix *P0;

	marray3d *alphaSamples;
	marray3d *Z;
	marray3d *H;
	marray3d *T;
	marray3d *R;
	marray3d *Q;

	_NGPsample *sampler;

} NGPsample;

typedef struct{
	int Nstates;
	int Nsamples;
	int approx;
	int Np;
	int Ns;

	gsl_vector *a0;

	gsl_matrix *P0;

	marray3d *Zt;
	marray3d *Ht;
	marray3d *R;
	marray3d *T;
	marray3d *Q;
	marray3d *alphaSamples;

	_NGPsample *sampler;

} NGPsampleZH;

/***************************************
* Constructor, destructor, main method *
***************************************/

_NGPsample * _NGPsample_New();
int _NGPsample_init(_NGPsample *self);
int _NGPsample_free(_NGPsample *s);
int _NGPsample_construct(_NGPsample *self,gsl_matrix *y,int Ns,gsl_vector *a0,gsl_matrix *P0,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,marray3d *alphaSamples);
int _NGPsample_operations(_NGPsample *self);

/****************************************
* Class method for doing operation once *
****************************************/

int _NGPsample_Sample(gsl_matrix *y,int Ns,gsl_vector *a0,gsl_matrix *P0,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,marray3d *alphaSamples);

/*********************************
* Writing method outputs to file *
*********************************/

int _NGPsample_writeAlphaSamples(_NGPsample *self,char *name);
int _NGPsample_writeOutputs(_NGPsample *self,char *baseName);


NGPsample * NGPsample_New();
int NGPsample_init(NGPsample *self);
int NGPsample_free(NGPsample *s);
int NGPsample_construct(NGPsample *self,gsl_matrix *y,int Nsamples,int Ns,gsl_vector *delta,double *sigEps,gsl_vector *sigU,gsl_vector *sigA,double *sigMu,double *sigAlph,int approx);
int NGPsample_operations(NGPsample *self);
//int NGPsample_update(NGPsample *self,gsl_matrix *y,int Nsamples,int Ns,gsl_vector *delta,double *sigEps,gsl_vector *sigU,gsl_vector *sigA,double *sigMu,double *sigAlph,int approx,marray *alphaSamples);

NGPsampleZH * NGPsampleZH_New();
int NGPsampleZH_init(NGPsampleZH *self);
int NGPsampleZH_free(NGPsampleZH *s);
int NGPsampleZH_construct(NGPsampleZH *self,gsl_matrix *y,int Nsamples,marray3d *Z,marray3d *H,gsl_vector *delta,gsl_vector *sigU,gsl_vector *sigA,double *sigMu,double *sigAlph,int approx,marray3d *alphaSamples);
int NGPsampleZH_operations(NGPsampleZH *self);


#endif
