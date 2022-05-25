#ifndef _ABT_OBJECTKALMAN_H
#define _ABT_OBJECTKALMAN_H

#include <gsl/gsl_matrix.h>
#include <marray.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
	gsl_matrix *HH;
	gsl_matrix *TT;
	gsl_matrix *RR;
	gsl_matrix *QQ;
	gsl_matrix *v;
	gsl_matrix *ZZ;
	gsl_matrix *y;
	gsl_matrix *ZPt_Mpm;
	gsl_matrix *Zt_Mmp;
	gsl_matrix *Pt_Mmm;
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
	marray *P_Atmm;
	marray *T;
	marray *R;
	marray *Q;
	marray *K;
	marray *Finv;
	marray *H;
	marray *Z;

} Kalman;

/***************************************
* Constructor, destroctor, main method *
***************************************/

Kalman * Kalman_New();
int Kalman_init(Kalman *self);
int Kalman_free(Kalman *s);
int Kalman_construct(Kalman *self,gsl_matrix *y,gsl_vector *a_init,gsl_matrix *P_init,marray *Z,marray *H,marray *T,marray *R,marray *Q,gsl_matrix *v,marray *K,marray *Finv);
int Kalman_operations(Kalman *self);

/****************************************
* Class method for doing operation once *
****************************************/

int Kalman_Filter(gsl_matrix *y,gsl_vector *a_init,gsl_matrix *P_init,marray *Z,marray *H,marray *T,marray *R,marray *Q,gsl_matrix *v,marray *K,marray *Finv);

/*********************************
* Writing method outputs to file *
*********************************/

int Kalman_writeV(Kalman *self,char *name);
int Kalman_writeK(Kalman *self,char *name);
int Kalman_writeFinv(Kalman *self,char *name);
int Kalman_writeOutputs(Kalman *self,char *baseName);

//Getters and setters of the model
//This comes exactly from Java
//Note there are no setters or getters for the outputs
//Also every time you set or get it zeros out whatever was
//in the output so we arent misled on the state
int Kalman_setY(Kalman *self,gsl_matrix *y);
int Kalman_getY(Kalman *self,gsl_matrix *y);

int Kalman_setA_init(Kalman *self,gsl_vector *a_init);
int Kalman_getA_init(Kalman *self,gsl_vector *a_init);

int Kalman_getP_init(Kalman *self,gsl_matrix *P0);
int Kalman_setP_init(Kalman *self,gsl_matrix *P0);

int Kalman_getZ(Kalman *self,marray *Z);
int Kalman_setZ(Kalman *self,marray *Z);

int Kalman_setH(Kalman *self,marray *H);
int Kalman_getH(Kalman *self,marray *H);

int Kalman_setT(Kalman *self,marray *T);
int Kalman_getT(Kalman *self,marray *T);

int Kalman_setR(Kalman *self,marray *R);
int Kalman_getR(Kalman *self,marray *R);

int Kalman_setQ(Kalman *self,marray *Q);
int Kalman_getQ(Kalman *self,marray *Q);

int Kalman_getV(Kalman *self,gsl_matrix *v);
int Kalman_getK(Kalman *self,marray *K);
int Kalman_getF(Kalman *self,marray *F);

#endif
