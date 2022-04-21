#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include "mdarray_complex.h"
#include "gmlib.h"
#include "KalmanFilter2c.h"
#include "matrix_utils.h"

static int _Kalman2c_checkSizes(Kalman2c *self);
static int _Kalman2c_zeroOutputs(Kalman2c *self);

Kalman2c * Kalman2c_New(){
	Kalman2c *self;
	self = (Kalman2c *)GM_Malloc(sizeof(Kalman2c));
	if(self==NULL)								GMERR(-11);
	if(Kalman2c_init(self))						GMERR(-21);
	return(self);
GMERRH("Kalman2c_New",NULL);
}

int Kalman2c_init(Kalman2c *self){
	if(self==NULL)								GMERR(-11);
	memset(self,'\0',sizeof(Kalman2c));

	return(0);
GMERRH("Kalman2c_init",1);
}

/*This frees all the temporary variables created. Note that
 *we DO NOT want to free the inputs here as that would cause
 *a problem if they were used by other objects.
*/ 

int Kalman2c_free(Kalman2c *s){
	if(s==NULL)									GMERR(-1);
	if(&GM_FreeGSLVector==NULL)					GMERR(-2);
	if(&GM_FreeGSLMatrix==NULL)					GMERR(-3);
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->yt_Vp,&s->at_Vm,&s->zat_Vp,&s->TTat_Vm,
				&s->KtVt_Vm,NULL))				GMERR(-21);
	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
				&s->a_Mtm,&s->RRQQ,
				&s->ZPt_Mpm,&s->Zt_Mmp,
				&s->ZPZ_Mpp,&s->Finvt_Mpp,&s->TTPt_Mmm,
				&s->TTPtZt_Mmp,&s->Kt_Mmp,&s->KtZ_Mmm,
				&s->L_Mmm,&s->Lt_Mmm,&s->TTPtLt_Mmm,
				&s->RQR,&s->RRT,NULL))			GMERR(-31);
	if(marray3d_complex_free(s->P_Atmm))		GMERR(-41);

	return(0);
GMERRH("Kalman2_free",1);
}

#define V_ALLOC(V,SIZE) self->V = gsl_vector_complex_alloc(SIZE)
#define M_ALLOC(M,SIZE1,SIZE2) self->M = gsl_matrix_complex_alloc(SIZE1,SIZE2)
#define A_ALLOC3(A,SIZE1,SIZE2,SIZE3) self->A = marray3d_complex_alloc(SIZE1,SIZE2,SIZE3)

int Kalman2c_construct(Kalman2c *self,gsl_matrix_complex *y,gsl_vector_complex *a_init,gsl_matrix_complex *P_init,marray3d_complex *Z,marray3d_complex *H,marray3d_complex *T,marray3d_complex *R,marray3d_complex *Q,gsl_matrix_complex *v,marray3d_complex *K,marray3d_complex *Finv){
	int t,p,m,r;
	self->Nt = y->size1;
    self->Np = y->size2;
    self->Nm = a_init->size;
	self->Nr = Q->d2;
	int H1 = H->d1;
	int T1 = T->d1;
	int R1 = R->d1;
	int Q1 = Q->d1;
	p = self->Np;
	m = self->Nm;
	r = self->Nr;
	t = self->Nt;

	//Allocate memory

	//Vectors
	self->a_init = a_init;
	V_ALLOC(yt_Vp,p);
	V_ALLOC(at_Vm,m);
	V_ALLOC(zat_Vp,p);
	V_ALLOC(TTat_Vm,m);
	V_ALLOC(KtVt_Vm,m);

	//Matrices
	self->P_init = P_init;
	self->v = v;
	self->y = y;
	M_ALLOC(a_Mtm,t,m);
	M_ALLOC(ZPt_Mpm,p,m);
	M_ALLOC(Zt_Mmp,m,p);
	M_ALLOC(ZPZ_Mpp,p,p);
	M_ALLOC(Finvt_Mpp,p,p);
	M_ALLOC(TTPt_Mmm,m,m);
	M_ALLOC(TTPtZt_Mmp,m,p);
	M_ALLOC(Kt_Mmp,m,p);
	M_ALLOC(KtZ_Mmm,m,m);
	M_ALLOC(L_Mmm,m,m);
	M_ALLOC(Lt_Mmm,m,m);
	M_ALLOC(TTPtLt_Mmm,m,m);
	M_ALLOC(RRQQ,m,r);
	M_ALLOC(RQR,m,m);
	M_ALLOC(RRT,r,m);


	//Arrays
	A_ALLOC3(P_Atmm,t,m,m);
	self->H = H;
	self->T = T;
	self->R = R;
	self->Q = Q;
	self->K = K;
	self->Finv = Finv;
	self->Z = Z;

	return(0);
GMERRH("Kalman2c_construct",1);
}

#undef M_ALLOC
#undef V_ALLOC
#undef A_ALLOC3

int Kalman2c_Filter(gsl_matrix_complex *y,gsl_vector_complex *a_init,gsl_matrix_complex *P_init,marray3d_complex *Z,marray3d_complex *H,marray3d_complex *T,marray3d_complex *R,marray3d_complex *Q,gsl_matrix_complex *v,marray3d_complex *K,marray3d_complex *Finv){
	Kalman2 *self = NULL;
	if((self=Kalman2c_New())==NULL)				GMERR(-11);
	if(Kalman2c_construct(self,y,a_init,P_init,Z,H,T,R,Q,v,K,Finv))	GMERR(-21);
	if(_Kalman2c_checkSizes(self))				GMERR(-30);
	if(Kalman2c_operations(self))				GMERR(-31);
	if(Kalman2c_free(self))						GMERR(-41);

	return(0);
GMERRH("Kalman2c_Filter",1);
}

/******************************************************
* Dimension checking for all the inputs and ouputs of *
* the kalman filter.								  *
******************************************************/

static int _Kalman2c_checkSizes(Kalman2c *self){
	int Np,Nm,Nr,Nt;
	if(self==NULL)											GMERR(-1);
	Np = self->Np;Nm = self->Nm;Nr = self->Nr;Nt = self->Nt;

	//Check inputs
	if(marray3d_complex_checkYZ(self->Z,Np,Nm))				GMERR(-11);
	if(marray3d_complex_checkYZ(self->H,Np,Np))				GMERR(-21);
	if(marray3d_complex_checkYZ(self->T,Nm,Nm))				GMERR(-31);
	if(marray3d_complex_checkYZ(self->R,Nm,Nr))				GMERR(-41);
	if(marray3d_complex_checkYZ(self->Q,Nr,Nr))				GMERR(-51);
	if(matrixCheckSize(self->v,Nt,Np))						GMERR(-61);
	if(matrixCheckSize(self->P_init,Nm,Nm))					GMERR(-71);
	if(marray3d_complex_checkX1orT(self->Z,Nt))				GMERR(-81);
	if(marray3d_complex_checkX1orT(self->H,Nt))				GMERR(-91);
	if(marray3d_complex_checkX1orT(self->T,Nt))				GMERR(-101);
	if(marray3d_complex_checkX1orT(self->R,Nt))				GMERR(-111);
	if(marray3d_complex_checkX1orT(self->Q,Nt))				GMERR(-121);

	//Check outputs
	if(marray3d_complex_checkYZ(self->Q,Nr,Nr))				GMERR(-51);
	if(marray3d_complex_checkYZ(self->Q,Nr,Nr))				GMERR(-51);
	
	return(0);
GMERRH("_Kalman2c_checkSizes",1);
}

/*******************************************************
* This method is called at the beginning of operations *
* to ensure that when something is changed the results *
* are from the current run. Redundant				   *
*******************************************************/

static int _Kalman2c_zeroOutputs(Kalman2c *self){
	if(self==NULL)											GMERR(-1);
	gsl_matrix_complex_set_zero(self->v);
	gsl_matrix_complex_set_zero(self->a_Mtm);
	if(marray3d_complex_set_zero(self->K))					GMERR(-11);
	if(marray3d_complex_set_zero(self->Finv))				GMERR(-21);
	if(marray3d_complex_set_zero(self->P_Atmm))				GMERR(-31);

	return(0);
GMERRH("_Kalman2c_zeroOutputs",1);
}
/*
#define A_GET_ROW(A,M,T) marray3d_get_X(self->A,self->M,T)
#define A_SET_ROW(A,M,T) marray3d_set_X(self->A,self->M,T)
#define DDOT_MV(M,V,R) gsl_blas_dgemv(CblasNoTrans,1.0,self->M,self->V,0.0,self->R)
#define DDOT_VV(UU,VV,RR) gsl_blas_ddot(self->UU,self->VV,self->RR)
#define DDOT_MM(MM,NN,RR) gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,self->MM,self->NN,0.0,self->RR)
#define M_GET_ROW(M,V,T) gsl_matrix_get_row(self->V,self->M,T)
#define M_SET_ROW(M,V,T) gsl_matrix_set_row(self->M,T,self->V)
#define M_ADD(MM,NN) gsl_matrix_add(self->MM,self->NN)
#define V_ADD(VV,UU) gsl_vector_add(self->VV,self->UU)
#define V_SUB(VV,UU) gsl_vector_sub(self->VV,self->UU)
#define M_SUB(MM,NN) gsl_matrix_sub(self->MM,self->NN)
#define A_ROW(A,M,T) M = self->A->matrixList[T]
*/

#define A_GET_ROW(A,M,T) marray3d_complex_get_X(self->A,self->M,T)
#define A_SET_ROW(A,M,T) marray3d_complex_set_X(self->A,self->M,T)
#define DDOT_MV(M,V,R) gsl_blas_zgemv(CblasNoTrans,1.0,self->M,self->V,0.0,self->R)
#define DDOT_VV(UU,VV,RR) gsl_blas_zdotu(self->UU,self->VV,self->RR)
#define DDOT_MM(MM,NN,RR) gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,1.0,self->MM,self->NN,0.0,self->RR)
#define M_GET_ROW(M,V,T) gsl_matrix_get_row(self->V,self->M,T)
#define M_SET_ROW(M,V,T) gsl_matrix_set_row(self->M,T,self->V)
#define M_ADD(MM,NN) gsl_matrix_add(self->MM,self->NN)
#define M_SUB(MM,NN) gsl_matrix_sub(self->MM,self->NN)
#define V_ADD(VV,UU) gsl_vector_add(self->VV,self->UU)
#define V_SUB(VV,UU) gsl_vector_sub(self->VV,self->UU)
#define A_ROW(A,M,T) M = self->A->matrixList[T]

/*
    Perform Kalman filtering on the data y.
    Conventions are as in Durbin and Koopman (2012).
    Relevant dimensions are:
    Nt: number of time points
    Np: dimension of observation space
    Nm: dimension of state space
    Nr: dimension of state noise covariance
    Parameters:
    a_init: mean of prior on states: (Nm, 1) or (Nm,)
    P_init: variance of prior on states: (Nm, Nm)
    Z: either a (Nt, Np, Nm) array or a (Np, Nm) array
    H: either a (Nt, Np, Np) array or a (Np, Np) array
    T: either a (Nt, Nm, Nm) array or a (Nm, Nm) array
    R: either a (Nt, Nm, Nr) array or a (Nm, Nr) array
    Q: either a (Nt, Nr, Nr) array or a (Nr, Nr) array
    Returns:
    v: (Nt, Np) vector of residuals at each time
    K: (Nt, Nm, Np) Kalman gain matrix
    Finv: (Nt, Np, Np) inverse of prediction variance matrix
*/

int Kalman2_operations(Kalman2 *self){
	int t;
	int Nt = self->Nt;
	//These are to see if we need to update state matrices every time
	int Ht = ((self->H->d1)>1);
	int Tt = ((self->T->d1)>1);
	int Rt = ((self->R->d1)>1);
	int Qt = ((self->Q->d1)>1);
	int Zt = ((self->Z->d1)>1);

	//Declare pointers to matrices
	gsl_matrix_complex *HH,*TT,*RR,*QQ,*ZZ,*Pt_Mmm;

	//Zero the outputs
	if(_Kalman2c_zeroOutputs(self))											GMERR(-11);
	
	//Initialize the variables
	if(M_SET_ROW(a_Mtm,a_init,0))											GMERR(-21);
	if(A_SET_ROW(P_Atmm,P_init,0))											GMERR(-31);
	A_ROW(H,HH,0);
	A_ROW(T,TT,0);
	A_ROW(R,RR,0);
	A_ROW(Q,QQ,0);
	A_ROW(Z,ZZ,0);
	A_ROW(P_Atmm,Pt_Mmm,0);
//	if(gTranspose(ZZ,self->Zt_Mmp))											GMERR(-41);
	gsl_matrix_complex_memcpy(self->Zt_Mmp,ZZ);
	
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,RR,QQ,0.0,self->RRQQ);
//	if(gTranspose(RR,self->RRT))											GMERR(-51);
	gsl_matrix_complex_memcpy(self->RRT,RR);
	for(t=0;t<Nt;t++){
		//Update the matrices if required

		if(Ht) A_ROW(H,HH,t);
		if(Tt) A_ROW(T,TT,t);
		if(Rt) A_ROW(R,RR,t);
		if(Qt) A_ROW(Q,QQ,t);
		if(Zt) A_ROW(Z,ZZ,t);
//		if(Zt) if(gTranspose(ZZ,self->Zt_Mmp))								GMERR(-132);
		if(Zt) gsl_matrix_complex_transpose_memcpy(self->Zt_Mmp,ZZ);
	
		//v[t] = y[t] - ZZ.dot(a[t])
		if(M_GET_ROW(a_Mtm,at_Vm,t))										GMERR(-141);
		if(M_GET_ROW(y,yt_Vp,t))											GMERR(-151);

		//if(DDOT_MV(ZZ,at_Vm,zat_Vp))										GMERR(-161);
		gsl_blas_dgemv(CblasNoTrans,1.0,ZZ,self->at_Vm,0.0,self->zat_Vp);
		if(V_SUB(yt_Vp,zat_Vp))												GMERR(-171);
		if(M_SET_ROW(v,yt_Vp,t))											GMERR(-181);

		//F = Z.dot(P[t]).dot(Z.T) + HH
		A_ROW(P_Atmm,Pt_Mmm,t);
		//if(DDOT_MM(ZZ,Pt_Mmm,ZPt_Mpm))										GMERR(-201);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,ZZ,Pt_Mmm,0.0,self->ZPt_Mpm);
		if(DDOT_MM(ZPt_Mpm,Zt_Mmp,ZPZ_Mpp))									GMERR(-211);
		//if(M_ADD(ZPZ_Mpp,HH))												GMERR(-221);
		if(gsl_matrix_add(self->ZPZ_Mpp,HH))								GMERR(-221);
		
		//Finv[t] = np.linalg.inv(F)
//		if(inv_LU(self->ZPZ_Mpp,self->Finvt_Mpp))							GMERR(-231);
		if(A_SET_ROW(Finv,Finvt_Mpp,t))										GMERR(-241);

		//K[t] = TT.dot(P[t]).dot(Z.T).dot(Finv[t])
		//if(DDOT_MM(TT,Pt_Mmm,TTPt_Mmm))										GMERR(-251);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,TT,Pt_Mmm,0.0,self->TTPt_Mmm);
		if(DDOT_MM(TTPt_Mmm,Zt_Mmp,TTPtZt_Mmp))								GMERR(-261);
		if(DDOT_MM(TTPtZt_Mmp,Finvt_Mpp,Kt_Mmp))							GMERR(-271);
		if(A_SET_ROW(K,Kt_Mmp,t))											GMERR(-281);

		//L = TT - K[t].dot(Z)
		//if(DDOT_MM(Kt_Mmp,ZZ,KtZ_Mmm))										GMERR(-291);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,self->Kt_Mmp,ZZ,0.0,self->KtZ_Mmm);
		if(gsl_matrix_memcpy(self->L_Mmm,TT))								GMERR(-301);
		if(M_SUB(L_Mmm,KtZ_Mmm))											GMERR(-311);

		if(t+1<Nt){
			//a[t+1] = TT.dot(a[t]) + K[t].dot(v[t])
			//if(DDOT_MV(TT,at_Vm,TTat_Vm))									GMERR(-321);
			gsl_blas_dgemv(CblasNoTrans,1.0,TT,self->at_Vm,0.0,self->TTat_Vm);

			if(DDOT_MV(Kt_Mmp,yt_Vp,KtVt_Vm))								GMERR(-331);
			if(V_ADD(TTat_Vm,KtVt_Vm))										GMERR(-341);
			if(M_SET_ROW(a_Mtm,TTat_Vm,t+1))								GMERR(-351);

			//P[t+1] = TT.dot(P[t]).dot(L.T) + RR.dot(QQ).dot(RR.T)
			if(gTranspose(self->L_Mmm,self->Lt_Mmm))						GMERR(-361);
			if(DDOT_MM(TTPt_Mmm,Lt_Mmm,TTPtLt_Mmm))							GMERR(-371);
	//		if(DDOT_MM(RR,QQ,RRQQ))											GMERR(-381);
			if(Rt||Qt) 	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,RR,QQ,0.0,self->RRQQ);
			if(Rt)		if(gTranspose(RR,self->RRT))						GMERR(-391);
			if(DDOT_MM(RRQQ,RRT,RQR))										GMERR(-401);
			if(M_ADD(TTPtLt_Mmm,RQR))										GMERR(-411);
			if(A_SET_ROW(P_Atmm,TTPtLt_Mmm,t+1))							GMERR(-421);
		}
	}
	return(0);
GMERRH("Kalman2c_operations",1);
}

#undef A_GET_ROW
#undef A_SET_ROW
#undef DDOT_MV
#undef DDOT_VV
#undef DDOT_MM
#undef M_GET_ROW
#undef M_SET_ROW
#undef M_ADD
#undef V_ADD
#undef M_SUB
#undef V_SUB
#undef A_ROW

/*********************************
* Writing method outputs to file *
*********************************/

int Kalman2c_writeV(Kalman2 *self,char *name){
	FILE *fp;
	if(self==NULL)					GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)					GMERR(-3);
	if(gsl_matrix_fprintf(fp,self->v,"%.15f"))	GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("Kalman2c_writeV",1);
}

int Kalman2c_writeK(Kalman2 *self,char *name){
	FILE *fp;
	if(self==NULL)					GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)					GMERR(-3);
	if(marray3d_write(fp,self->K))				GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("Kalman2c_writeK",1);
}

int Kalman2c_writeFinv(Kalman2 *self,char *name){
	FILE *fp;
	if(self==NULL)					GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)					GMERR(-3);
	if(marray3d_write(fp,self->Finv))			GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("Kalman2c_writeFinv",1);
}

int Kalman2c_writeOutputs(Kalman2 *self,char *baseName){
	//Here we make all the output names
	char strV[16] = "_V.txt";
	char strK[16] = "_K.txt";
	char strFinv[16] = "_Finv.txt";
	char outV[1024],outK[1024],outFinv[1024];
	sprintf(outV,"%s%s",baseName,strV);
	sprintf(outK,"%s%s",baseName,strK);
	sprintf(outFinv,"%s%s",baseName,strFinv);

	//Here we actually write out all the data
	if(self==NULL)						GMERR(-1);
	if(Kalman2c_writeV(self,outV))		GMERR(-11);
	if(Kalman2c_writeK(self,outK))		GMERR(-21);
	if(Kalman2c_writeFinv(self,outFinv))	GMERR(-31);

	return(0);
GMERRH("Kalman2c_writeOutputs",1);
}

