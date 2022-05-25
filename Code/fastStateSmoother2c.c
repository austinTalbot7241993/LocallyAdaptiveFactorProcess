#include <stdio.h>
#include <string.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix_complex.h>
#include <gsl/gsl_vector_complex.h>
#include <math.h>
#include "gmlib.h"
#include "matrix_utils.h"
#include "mdarray_complex.h"
#include "fastStateSmoother2c.h"

/********************************
*
*
*
********************************/

static int _FSS2_checkSizes(FSS2 *self);

FSS2 * FSS2_New(){
	FSS2 *self;
	self = (FSS2 *)GM_Malloc(sizeof(FSS2));
	if(self==NULL)							GMERR(-11);
	if(FSS2_init(self))						GMERR(-21);
	return(self);
GMERRH("FSS2_New",NULL);
}

int FSS2_init(FSS2 *self){
	if(self==NULL)							GMERR(-11);
	memset(self,'\0',sizeof(FSS2));

	return(0);
GMERRH("FSS2_init",1);
}


int FSS2_free(FSS2 *s){
	if(s==NULL)							GMERR(-11);
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->vt_Vp,&s->Fvt_Vp,&s->Kttrt_Vp,
				&s->rt_Vm,&s->r_init,&s->TTrt_Vm,
				&s->ZZFvt_Vm,&s->Pr_Vm,&s->TTat_Vm,
				&s->RQRrt_Vm,&s->at_Vm,NULL))	GMERR(-21);
	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
				&s->r,&s->ZZT_Mmp,&s->RRQQ,&s->RRt,
				&s->RQR,&s->Ktt_Mpm,
				&s->TTt_Mmm,NULL))	GMERR(-31);

	return(0);
GMERRH("FSS2_free",1);
}

#define M_ALLOC(M,SIZE1,SIZE2) self->M = gsl_matrix_alloc(SIZE1,SIZE2)
#define V_ALLOC(V,SIZE) self->V = gsl_vector_alloc(SIZE)

int FSS2_construct(FSS2 *self,gsl_matrix *v,marray3d *K,marray3d *Finv,gsl_vector *a_init,gsl_matrix *P_init,marray3d *Z,marray3d *T,marray3d *R,marray3d *Q,gsl_matrix *alpha){
	int t,p,m,r;
	self->Nt = v->size1;
	self->Np = v->size2;
	self->Nm = a_init->size;
	//self->Nr = getSizeY(Q);
	self->Nr = Q->d2;

	t = self->Nt;
	p = self->Np;
	m = self->Nm;
	r = self->Nr;

	//Allocate memory

	//Vectors
	//Input
	self->a_init = a_init;

	//Temporary
	V_ALLOC(vt_Vp,p);
	V_ALLOC(Fvt_Vp,p);
	V_ALLOC(Kttrt_Vp,p);
	V_ALLOC(rt_Vm,m);
	V_ALLOC(r_init,m);
	V_ALLOC(TTrt_Vm,m);
	V_ALLOC(ZZFvt_Vm,m);
	V_ALLOC(Pr_Vm,m);
	V_ALLOC(TTat_Vm,m);
	V_ALLOC(RQRrt_Vm,m);
	V_ALLOC(at_Vm,m);


	//Matices
	//Input
	self->v = v;
	self->P_init = P_init;
	self->alpha=alpha;

	//Temporary
	M_ALLOC(r,t,m);
	M_ALLOC(ZZT_Mmp,m,p);
	M_ALLOC(RRQQ,m,r);
	M_ALLOC(RRt,r,m);
	M_ALLOC(RQR,m,m);
	M_ALLOC(Ktt_Mpm,p,m);
	M_ALLOC(TTt_Mmm,m,m);

	//Arrays
	//Input
	self->Z = Z;
	self->K = K;
	self->Finv = Finv;
	self->T = T;
	self->R = R;
	self->Q = Q;

	//Temporary
	return(0);
GMERRH("FSS2_construct",1);
}

#undef M_ALLOC
#undef V_ALLOC

/*This method ensures that all the inputs are feasible sizes to actually
 *create a state smoother */

static int _FSS2c_checkSizes(FSS2c *self){
	int Np,Nm,Nr,Nt;
	if(self==NULL)										GMERR(-1);
	Np=self->Np;Nm=self->Nm;Nr=self->Nr;Nt=self->Nt;

	if(marray3d_complex_checkSizes(self->K,Nt,Nm,Np))	GMERR(-11);
	if(marray3d_complex_checkSizes(self->Finv,Nt,Np,Np))GMERR(-21);
	if(marray3d_complex_checkYZ(self->T,Nm,Nm))			GMERR(-31);
	if(marray3d_complex_checkYZ(self->Q,Nr,Nr))			GMERR(-41);
	if(marray3d_complex_checkYZ(self->Z,Np,Nm))			GMERR(-51);

	if(matrixComplexCheckSize(self->v,Nt,Np))			GMERR(-61);
	if(matrixComplexCheckSize(self->P_init,Nm,Nm))		GMERR(-71);
	if(matrixComplexCheckSize(self->alpha,Nt,Nm))		GMERR(-81);

	if(self->a_init->size!=Nm)							GMERR(-91);

	if(marray3d_complex_checkX1orT(self->T,Nt))			GMERR(-101);
	if(marray3d_complex_checkX1orT(self->R,Nt))			GMERR(-111);
	if(marray3d_complex_checkX1orT(self->Q,Nt))			GMERR(-121);
	if(marray3d_complex_checkX1orT(self->Z,Nt))			GMERR(-131);
	
	return(0);
GMERRH("_FSS2c_checkSizes",1);
}

/*This zeros out the quantities of interest in our methods so that our results
 *definitely aren't contaminated by a previous run.
*/

static int _FSS2_zeroOutputs(FSS2 *self){
	gsl_matrix_complex_set_zero(self->alpha);
	gsl_matrix_complex_set_zero(self->r);
	gsl_vector_complex_set_zero(self->r_init);
	return(0);
GMERRH("_FSS2c_zeroOutputs",1);
}



/*This is the class method for smoothing the inputs. Notice that if we are going to do this 
 *repeatedly you will want to create the object separately so that you don't have to repeat
 *the initialization and free steps*/

int FSS2_Smoother(gsl_matrix *v,marray3d *K,marray3d *Finv,gsl_vector *a_init,gsl_matrix *P_init,marray3d *Z,marray3d *T,marray3d *R,marray3d *Q,gsl_matrix *alpha){
	FSS2 *self = NULL;
	if((self=FSS2_New())==NULL)										GMERR(-11);
	if(FSS2_construct(self,v,K,Finv,a_init,P_init,Z,T,R,Q,alpha))	GMERR(-21);
	if(_FSS2_checkSizes(self))										GMERR(-31);
	if(FSS2_operations(self))										GMERR(-41);
	if(FSS2_free(self))												GMERR(-51);
	return(0);
GMERRH("FSS2_Smoother",1);

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
#define V_SUB(VV,UU) gsl_vector_sub(self->VV,self->UU)
#define V_ADD(VV,UU) gsl_vector_add(self->VV,self->UU)
#define M_SUB(MM,NN) gsl_matrix_sub(self->MM,self->NN)
#define A_ROW(A,M,T) M = self->A->matrixList[T]*/

#define A_GET_ROW(A,M,T) marray3d_complex_get_X(self->A,self->M,T)
#define A_SET_ROW(A,M,T) marray3d_complex_get_X(self->A,self->M,T)
#define DDOT_MV(M,V,R) gsl_blas_zgemv(CblasNoTrans,1.0,self->m,self->V,0.0,self->R)
#define DDOT_VV(UU,VV,RR) gsl_blas_zdotu(self->UU,self->VV,self->RR)
#define DDOT_MM(MM,NN,RR) gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,1.0,self->MM,self->NN,0.0,self->RR)
#define M_GET_ROW(M,V,T) gsl_matrix_complex_get_row(self->V,self->M,T)
#define M_SET_ROW(M,V,T) gsl_matrix_complex_set_row(self->M,T,self->V)
#define M_ADD(MM,NN) gsl_matrix_complex_add(self->MM,self->NN)
#define V_SUB(VV,UU) gsl_vector_complex_sub(self->VV,self->UU)
#define V_ADD(VV,UU) gsl_vector_complex_add(self->VV,self->UU)
#define M_SUB(MM,NN) gsl_matrix_complex_sub(self->MM,self->NN)
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

int FSS2_operations(FSS2 *self){
	int t;
	int Nt = self->Nt;
	int Tt = ((self->T->d1)>1);
	int Rt = ((self->R->d1)>1);
	int Qt = ((self->Q->d1)>1);
	int Zt = ((self->Z->d1)>1);

	gsl_matrix *ZZ,*TT,*RR,*QQ,*Kt_Mmp,*Finvt_Mpp;

	//Make sure that the outputs are uncontaminated by previous runs
	if(_FSS2_zeroOutputs(self))										GMERR(-1);

	//If single dimension don't update
	//if(A_GET_ROW(Z,ZZ,0))											GMERR(-11);
	A_ROW(Z,ZZ,0);
	A_ROW(T,TT,0);
	A_ROW(R,RR,0);
	A_ROW(Q,QQ,0);
	//if(gTranspose(ZZ,self->ZZT_Mmp))								GMERR(-12);
	gsl_matrix_complex_transpose_memcpy(self->ZZT_Mmp,ZZ);

	for(t=Nt-1;t>-1;t--){
		//Update the states for each time point if neccessary
		if(Tt) A_ROW(T,TT,t);
		if(Rt) A_ROW(R,RR,t);
		if(Qt) A_ROW(Q,QQ,t);
		if(Zt) A_ROW(Z,ZZ,t);
		if(Zt)	gsl_matrix_complex_transpose_memcpy(self->ZZT_Mmp,ZZ);

		//uu=Finv[t].dot(v[t])-K[t].T.dot(r[t])
		A_ROW(Finv,Finvt_Mpp,t);
		A_ROW(K,Kt_Mmp,t);
		if(M_GET_ROW(r,rt_Vm,t))									GMERR(-101);
		if(M_GET_ROW(v,vt_Vp,t))									GMERR(-111);
//		if(DDOT_MV(Finvt_Mpp,vt_Vp,Fvt_Vp))							GMERR(-121);
		if(gsl_blas_dgemv(CblasNoTrans,1.0,Finvt_Mpp,self->vt_Vp,0.0,self->Fvt_Vp))	GMERR(-121);
//		if(gTranspose(Kt_Mmp,self->Ktt_Mpm))						GMERR(-131);
		gsl_matrix_complex_transpose_memcpy(self->Ktt_Mpm,self->Kt_Mmp);
		if(DDOT_MV(Ktt_Mpm,rt_Vm,Kttrt_Vp))							GMERR(-141);
		if(V_SUB(Fvt_Vp,Kttrt_Vp))									GMERR(-151);

		//thisr = ZZ.T.dot(u) + TT.T.dot(r[t])
		if(DDOT_MV(ZZT_Mmp,Fvt_Vp,ZZFvt_Vm))						GMERR(-161);
//		if(gTranspose(self->TT,self->TTt_Mmm))						GMERR(-171);
//		if(gTranspose(TT,self->TTt_Mmm))							GMERR(-171);
		gsl_matrix_complex_transpose_memcpy(self->TTt_Mmm,TT);

		if(DDOT_MV(TTt_Mmm,rt_Vm,TTrt_Vm))							GMERR(-181);
		if(V_ADD(ZZFvt_Vm,TTrt_Vm))									GMERR(-191);
		if(t>0){
			//r[t-1] = thisr
			if(M_SET_ROW(r,ZZFvt_Vm,t-1))							GMERR(-201);
		}else{
			//r_init = thisr
			if(gsl_vector_memcpy(self->r_init,self->ZZFvt_Vm))		GMERR(-211);
		}
	}

	//alpha[0] = a_init + P_init.dot(r_init)
	if(DDOT_MV(P_init,r_init,Pr_Vm))									GMERR(-221);
	if(V_ADD(Pr_Vm,a_init))											GMERR(-231);
	if(M_SET_ROW(alpha,Pr_Vm,0))									GMERR(-241);

	//RQR = RR.dot(QQ).dot(RR.T)
	//if(DDOT_MM(RR,QQ,RRQQ))											GMERR(-251);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,RR,QQ,0.0,self->RRQQ))	GMERR(-251);
//	if(gTranspose(RR,self->RRt))									GMERR(-261);
	gsl_matrix_complex_transpose_memcpy(self->RRt,self->RR);
	if(DDOT_MM(RRQQ,RRt,RQR))										GMERR(-271);
	
	for(t=0;t<Nt-1;t++){
		//Update QQ,RR, RQR if required
		if(Tt) A_ROW(T,TT,t);
		if(Rt) A_ROW(R,RR,t);
		if(Qt) A_ROW(Q,QQ,t);
		
//		if(Rt||Qt) if(DDOT_MM(RR,QQ,RRQQ))							GMERR(-451);
		if(Rt||Qt) if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,RR,QQ,0.0,self->RRQQ))GMERR(-451);
		//if(Rt||Qt) if(gTranspose(RR,self->RRt))						GMERR(-461);
		if(Rt||Qt) gsl_matrix_complex_transpose_memcpy(self->RRt,self->RR);
		if(Rt||Qt) if(DDOT_MM(RRQQ,RRt,RQR))						GMERR(-471);

		//alpha[t+1] = TT.dot(alpha[t]) + RQR.dot(r[t])
		if(M_GET_ROW(alpha,at_Vm,t))								GMERR(-281);
		if(M_GET_ROW(r,rt_Vm,t))									GMERR(-291);
//		if(DDOT_MV(TT,at_Vm,TTat_Vm))								GMERR(-301);
		if(gsl_blas_dgemv(CblasNoTrans,1.0,TT,self->at_Vm,0.0,self->TTat_Vm)) GMERR(-121)
		if(DDOT_MV(RQR,rt_Vm,RQRrt_Vm))								GMERR(-311);
		if(V_ADD(TTat_Vm,RQRrt_Vm))									GMERR(-321);
		if(M_SET_ROW(alpha,TTat_Vm,(t+1)))							GMERR(-331);
	}
	return(0);
GMERRH("FSS2_operations",1);
}

#undef A_GET_ROW
#undef A_SET_ROW
#undef M_GET_ROW
#undef M_SET_ROW
#undef M_ADD
#undef M_SUB
#undef V_ADD
#undef V_SUB
#undef DDOT_MV
#undef DDOT_MM
#undef DDOT_VV
#undef A_ROW


int FSS2c_writeAlpha(FSS2 *self,char *name){
	FILE *fp;
	if(self==NULL)				GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)				GMERR(-3);
	if(gsl_matrix_complex_fprintf(fp,self->alpha,"%.15f"))			GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("FSS2c_writeAlpha",1);
}

int FSS2c_writeOutputs(FSS2 *self,char *baseName){
	FILE *fp;
	char strA[16] = "_Alpha.txt";
	char outA[1024];
	sprintf(outA,"%s%s",baseName,strA);
	if(self==NULL)					GMERR(-1);
	if(FSS2c_writeAlpha(self,outA))	GMERR(-11);
	return(0);
GMERRH("FSS2c_writeOutputs",1);
}

