#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "matrix_utils.h"
#include "gmlib.h"
#include "KalmanFilter2.h"
#include "fastStateSmoother2.h"
#include "SSsimulate2.h"
#include "NGPtools2.h"
#include "NGPsample2.h"
#include "mdarray.h"
#include "matrix_utils.h"
#include "mvnorm.h"

static int _NGPsample_checkSizes(_NGPsample *self,gsl_matrix *y,int Ns,gsl_vector *a0,gsl_matrix *P0,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,marray3d *alphaSamples);

_NGPsample * _NGPsample_New(){
	_NGPsample *self;
	self = (_NGPsample *)GM_Malloc(sizeof(_NGPsample));
	if(self==NULL)							GMERR(-11);
	if(_NGPsample_init(self))				GMERR(-21);
	return(self);
GMERRH("_NGPsample2_New",NULL);
}

int _NGPsample_init(_NGPsample *self){
	if(self==NULL)							GMERR(-11);
	memset(self,'\0',sizeof(_NGPsample));

	return(0);
GMERRH("_NGPsample2_init",1);
}

#undef NNN

int _NGPsample_free(_NGPsample *s){
	if(s==NULL)								GMERR(-1);
	if(SSsimulate2_free(s->SSsample))		GMERR(-11);
	gsl_matrix_free(s->alphaDraw);

	return(0);
GMERRH("_NGPsample2_free",1);
}

static int _NGPsample_checkSizes(_NGPsample *self,gsl_matrix *y,int Ns,gsl_vector *a0,gsl_matrix *P0,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,marray3d *alphaSamples){
	int m,p,r,s,t;
	if(self==NULL)							GMERR(-1);
	m = self->Nm;
	p = self->Np;
	t = self->Nt;
	r = self->Nr;
	s = self->Ns;

	if(marray3d_checkSizes(alphaSamples,s,t,m))						GMERR(-11);
	if(a0->size!=m)													GMERR(-21);
	if(matrixCheckSize(P0,m,m))										GMERR(-31);
	if(marray3d_checkYZ(Z,p,m))										GMERR(-41);
	if(marray3d_checkYZ(H,p,p))										GMERR(-51);
	if(marray3d_checkYZ(T,m,m))										GMERR(-61);
	if(marray3d_checkYZ(R,m,r))										GMERR(-71);
	if(marray3d_checkYZ(Q,r,r))										GMERR(-81);
	if(matrixCheckSize(y,t,p))										GMERR(-91);

	if(marray3d_checkX1orT(Z,t))									GMERR(-101);
	if(marray3d_checkX1orT(H,t))									GMERR(-111);
	if(marray3d_checkX1orT(T,t))									GMERR(-121);
	if(marray3d_checkX1orT(R,t))									GMERR(-131);
	if(marray3d_checkX1orT(Q,t))									GMERR(-141);
	
	return(0);
GMERRH("_NGPsample2_checkSizes",1);
}

//function _sample(y, Ns, a0, P0, Z, H, T, R, Q; interleaved=true)
#define M_ALLOC(M,SIZE1,SIZE2) self->M = gsl_matrix_alloc(SIZE1,SIZE2)

int _NGPsample_construct(_NGPsample *self,gsl_matrix *y,int Ns,gsl_vector *a0,gsl_matrix *P0,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,marray3d *alphaSamples){
	if(self==NULL)							GMERR(-1);
	self->Ns = Ns;
	self->Nr = Q->d2;//getSizeY(Q);
	self->Nt = y->size1;
	self->Np = y->size2;
	self->Nm = Z->d3;//getSizeZ(Z);
	if(_NGPsample_checkSizes(self,y,Ns,a0,P0,Z,H,T,R,Q,alphaSamples))	GMERR(-11);

	//Vectors
	self->a0 = a0;

	//Matrices
	self->P0 = P0;
	M_ALLOC(alphaDraw,self->Nt,self->Nm);

	//Marrays 
	self->Z = Z;
	self->H = H;
	self->T = T;
	self->R = R;
	self->Q = Q;
	self->alphaSamples = alphaSamples;

	
	self->SSsample = SSsimulate2_New();
	if(SSsimulate2_construct(self->SSsample,y,a0,P0,Z,H,T,R,Q,self->alphaDraw))	GMERR(-11);
	return(0);
GMERRH("_NGPsample2_construct",1);
}

#undef M_ALLOC

static int __NGPsample_checkSizes(_NGPsample *self){
	
	return(0);
GMERRH("__NGPsample2_checkSizes",1);
}

static int __NGPsample_zeroOutputs(_NGPsample *self){
	if(self==NULL)								GMERR(-1);
	if(self->alphaSamples==NULL)				GMERR(-2);
	if(self->alphaDraw==NULL)					GMERR(-3);
	marray3d_set_zero(self->alphaSamples);
	gsl_matrix_set_zero(self->alphaDraw);

	return(0);
GMERRH("__NGPsample2_zeroOutputs",1);
}

#define A_GET_ROW(A,M,T) marray3d_get_X(self->A,self->M,T)
#define A_SET_ROW(A,M,T) marray3d_set_X(self->A,self->M,T)

int _NGPsample_operations(_NGPsample *self){
	int idx,Ns;
	if(self==NULL)							GMERR(-1);
	if(__NGPsample_zeroOutputs(self))		GMERR(-11);
	Ns = self->Ns;
	for(idx=0;idx<Ns;idx++){
		//α_samples[:, :, idx] = simulate(y, a0, P0, Z, H, T, R, Q, interleaved=interleaved)
		if(SSsimulate2_operations(self->SSsample))							GMERR(-21);
		if(A_SET_ROW(alphaSamples,SSsample->alpha_draw,idx))				GMERR(-31);
	}

	return(0);
GMERRH("_NGPsample2_operations",1);
}

#undef A_GET_ROW
#undef A_SET_ROW

int _NGPsample_Sample(gsl_matrix *y,int Ns,gsl_vector *a0,gsl_matrix *P0,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,marray3d *alphaSamples){
	_NGPsample *self = NULL;
	if((self=_NGPsample_New())==NULL)		GMERR(-11);
	if(_NGPsample_construct(self,y,Ns,a0,P0,Z,H,T,R,Q,alphaSamples))	GMERR(-21);
	if(__NGPsample_checkSizes(self))		GMERR(-31);
	if(_NGPsample_operations(self))			GMERR(-41);
	if(_NGPsample_free(self))				GMERR(-51);
	
	return(0);
GMERRH("_NGPsample2_Sample",1);
}

/*********************************
* Writing method outputs to file *
*********************************/

int _NGPsample_writeAlphaSamples(_NGPsample *self,char *name){
	FILE *fp;
	if(self==NULL)							GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)							GMERR(-3);
	//if(marray_fprintf(fp,self->alphaSamples,"%.15f"))	GMERR(-11);
	if(marray3d_write(fp,self->alphaSamples))	GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("_NGPsample_writeAlphaSamples",1);
}
int _NGPsample_writeOutputs(_NGPsample *self,char *baseName){
	char strA[32] = "_AlphaSamples.txt";
	char outA[1024];
	sprintf(outA,"%s%s",baseName,strA);
	if(self==NULL)							GMERR(-1);
	if(_NGPsample_writeAlphaSamples(self,outA))	GMERR(-11);
	return(0);
GMERRH("_NGPsample_writeOutputs",1);
}

NGPsample * NGPsample_New(){
	NGPsample *self;
	self = (NGPsample *)GM_Malloc(sizeof(NGPsample));
	if(self==NULL)							GMERR(-11);
	if(NGPsample_init(self))				GMERR(-21);
	return(self);
GMERRH("NGPsample_New",NULL);
}

#define NNN(A) self->A=NULL
int NGPsample_init(NGPsample *self){
	if(self==NULL)							GMERR(-11);
	memset(self,'\0',sizeof(NGPsample));

	return(0);
GMERRH("NGPsample_init",1);
}
#undef NNN

int NGPsample_free(NGPsample *s){
	if(s==NULL)											GMERR(-1);
	if(_NGPsample_free(s->sampler))						GMERR(-11);
	gsl_vector_free(s->a0);
	gsl_matrix_free(s->P0);

	if(GM_FreezeListWithMethod(&GM_FreeGSLMarray,"x",
				&s->Z,&s->H,&s->T,&s->Q,&s->R,NULL))	GMERR(-21);
	
	return(0);
GMERRH("NGPsample_free",1);
}

#define M_ALLOC(M,SIZE1,SIZE2) self->M = gsl_matrix_alloc(SIZE1,SIZE2)
#define V_ALLOC(V,SIZE) self->V = gsl_vector_alloc(SIZE);
#define A_ALLOC3(A,SIZE1,SIZE2,SIZE3) self-> A = marray3d_alloc(SIZE1,SIZE2,SIZE3)

int NGPsample_construct(NGPsample *self,gsl_matrix *y,int Nsamples,int Ns,gsl_vector *delta,double *sigEps,gsl_vector *sigU,gsl_vector *sigA,double *sigMu,double *sigAlph,int approx){
	int Nm,Nr,Nt,Np;
	if(self==NULL)							GMERR(-1);
	Np = y->size2;
	self->Np = y->size2;
	self->Ns = Ns;
	Nt = y->size1;
	if(approx){
		Nr = 2 * Ns;
	}else{
		Nr = 3 * Ns;
	}
	Nm = 3*Ns;

	A_ALLOC3(Z,1,Np,3*Ns);
	A_ALLOC3(H,1,Ns,Ns);
	A_ALLOC3(T,Nt,Nm,Nm);
	A_ALLOC3(Q,Nt,Nr,Nr);
	A_ALLOC3(R,1,Nm,Nr);
	V_ALLOC(a0,Nm);
	M_ALLOC(P0,Nm,Nm);

	
	self->sampler = _NGPsample_New();
	if(NGP_assemble_matrices(self->Np,self->Ns,delta,sigEps,sigU,sigA,sigMu,sigAlph,approx,self->Z,self->H,self->T,self->R,self->Q,self->a0,self->P0))	GMERR(-11);
	if(_NGPsample_construct(self->sampler,y,self->Ns,self->a0,self->P0,self->Z,self->H,self->T,self->R,self->Q,self->alphaSamples))		GMERR(-21);

	return(0);
GMERRH("NGPsample_construct",1);
}

#undef M_ALLOC
#undef V_ALLOC
#undef A_ALLOC3
/*
int NGPsample_update(NGPsample *self,gsl_matrix *y,int Nsamples,int Ns,gsl_vector *delta,double *sigEps,gsl_vector *sigU,gsl_vector *sigA,double *sigMu,double *sigAlph,int approx){
	if(self==NULL)									GMERR(-1);
	if(NGP_assemble_matrices(self->Np,self->Ns,delta,sigEps,sigU,sigA,sigMu,sigAlph,approx,self->Z,self->H,self->T,self->R,self->Q,self->a0,self->P0))	GMERR(-11);
	if(_NGPsample_construct(self->sampler,y,self->Ns,self->a0,self->P0,self->Z,self->H,self->T,self->R,self->Q,self->alphaSamples))		GMERR(-21);

	return(0);
GMERRH("NGPsample_modify",1);
}
*/
int NGPsample_operations(NGPsample *self){
	if(self==NULL)			GMERR(-1);
	if(_NGPsample_operations(self->sampler))	GMERR(-11);
	return(0);
GMERRH("NGPsample_operations",1);
}

NGPsampleZH * NGPsampleZH_New(){
	NGPsampleZH *self;
	self = (NGPsampleZH *)GM_Malloc(sizeof(NGPsampleZH));
	if(self==NULL)							GMERR(-11);
	if(NGPsampleZH_init(self))				GMERR(-21);
	return(self);
GMERRH("NGPsampleZH_New",NULL);
}

int NGPsampleZH_init(NGPsampleZH *self){
	if(self==NULL)							GMERR(-11);
	memset(self,'\0',sizeof(NGPsampleZH));

	return(0);
GMERRH("NGPsampleZH_init",1);
}

int NGPsampleZH_free(NGPsampleZH *s){
	if(s==NULL)							GMERR(-1);
	if(GM_FreezeListWithMethod(&GM_FreeGSLMarray,"x",
				&s->T,&s->Q,&s->R,&s->Zt,&s->Ht,NULL))	GMERR(-21);
	gsl_vector_free(s->a0);
	gsl_matrix_free(s->P0);
	if(_NGPsample_free(s->sampler))				GMERR(-11);

	return(0);
GMERRH("NGPsampleZH_free",1);
}

int NGPsampleZH_construct(NGPsampleZH *self,gsl_matrix *y,int Nsamples,marray3d *Z,marray3d *H,gsl_vector *delta,gsl_vector *sigU,gsl_vector *sigA,double *sigMu,double *sigAlph,int approx,marray3d *alphaSamples){
	double eps = 0.0;
	double *sigEps = &eps;
	int Nstates = Z->d3;//getSizeZ(Z);
	if(self==NULL)							GMERR(-1);
	self->Nsamples = Nsamples;
	self->Np = y->size2;
	self->approx = approx;
	self->sampler = _NGPsample_New();

	if(NGP_assemble_matrices(self->Np,self->Ns,delta,sigEps,sigU,sigA,sigMu,sigAlph,approx,self->Zt,self->Ht,self->T,self->R,self->Q,self->a0,self->P0))  GMERR(-11);
	if(_NGPsample_construct(self->sampler,y,Nsamples,self->a0,self->P0,Z,H,self->T,self->R,self->Q,self->alphaSamples))	GMERR(-21);
	
	return(0);
GMERRH("NGPsampleZH_construct",1);
}

int NGPsampleZH_operations(NGPsampleZH *self){
	if(self==NULL)							GMERR(-1);
    if(_NGPsample_operations(self->sampler))    GMERR(-11);

	return(0);
GMERRH("NGPsampleZH_operations",1);
}


int NGPsampleZH_modify(NGPsampleZH *self,gsl_matrix *y,int Nsamples,marray3d *Z,marray3d *H,gsl_vector *delta,gsl_vector *sigU,gsl_vector *sigA,double *sigMu,double *sigAlph,int approx,marray3d *alphaSamples){
	double eps = 0.0;
	double *sigEps = &eps;
	int Nstates = Z->d3;//getSizeZ(Z);
	if(self==NULL)											GMERR(-1);

	if(NGP_assemble_matrices(self->Np,self->Ns,delta,sigEps,sigU,sigA,sigMu,sigAlph,approx,self->Zt,self->Ht,self->T,self->R,self->Q,self->a0,self->P0))  GMERR(-11);
	if(_NGPsample_construct(self->sampler,y,Nsamples,self->a0,self->P0,Z,H,self->T,self->R,self->Q,self->alphaSamples))	GMERR(-21);
	
	return(0);
GMERRH("NGPsampleZH_modify",1);
}
