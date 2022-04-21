#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "gmlib.h"
#include "matrix_utils.h"
#include "mdarray.h"
#include "KalmanFilter2.h"
#include "fastStateSmoother2.h"
#include "SSsimulate2.h"
#include <pthread.h>

static int _SSsimulate2_checkSizes(SSsimulate2 *self);
static int _SSsimulate2_zeroOutputs(SSsimulate2 *self);
int SSsimulate2_runParallel(SSsimulate2 *self);
static void *_run1(void *_self);
static void *_run2(void *_self);

/*****************************
******************************
**							**
** KFS substructure methods **
**							**
******************************
*****************************/

KFS * KFS_New(){
	KFS *self;
	self = (KFS *)GM_Malloc(sizeof(KFS));
	if(self==NULL)										GMERR(-11);
	if(KFS_init(self))									GMERR(-21);
	return(self);
GMERRH("KFS_New",NULL);
}


int KFS_init(KFS *self){
	if(self==NULL)										GMERR(-11);
	memset(self,'\0',sizeof(KFS));

	return(0);
GMERRH("KFS_init",1);
}

int KFS_free(KFS *s){
	if(s==NULL)											GMERR(-1);
	if(s->fss1!=NULL){
		if(FSS2_free(s->fss1))							GMERR(-11);
		s->fss1 = NULL;
	}
	if(s->kalman1!=NULL){
		if(Kalman2_free(s->kalman1))					GMERR(-21);
		s->kalman1 = NULL;
	}
	GM_Free(s);

	return(0);
GMERRH("KFS_free",1);
}

int KFS_construct(KFS *self,gsl_matrix *y,gsl_vector *a_init,gsl_matrix *P_init,
					marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,
					gsl_matrix *v,marray3d *K,marray3d *Finv1,gsl_matrix *alpha_hat){
	if(self==NULL)										GMERR(-1);
	self->kalman1 = Kalman2_New();
	self->fss1 = FSS2_New();
	printf(">>>>>10\n");
	if(Kalman2_construct(self->kalman1,y,a_init,P_init,Z,H,T,R,Q,v,K,Finv1))	GMERR(-11);
	printf(">>>>>20\n");
	if(FSS2_construct(self->fss1,v,K,Finv1,a_init,P_init,Z,T,R,Q,alpha_hat))	GMERR(-21);
	printf(">>>>>30\n");
	return(0);
GMERRH("KFS_construct",1);
}

int KFS_operations(KFS *self){
	if(self==NULL)										GMERR(-1);
		if(Kalman2_operations(self->kalman1))			GMERR(-11);
		if(FSS2_operations(self->fss1))					GMERR(-21);
	return(0);
GMERRH("KFS_operations",1);
}

/*
**
** SS simulation methods **
**
**/



SSsimulate2 * SSsimulate2_New(){
	SSsimulate2 *self;
	self = (SSsimulate2 *)GM_Malloc(sizeof(SSsimulate2));
	if(self==NULL)										GMERR(-11);
	if(SSsimulate2_init(self))							GMERR(-21);
	return(self);
GMERRH("SSsimulate2_New",NULL);
}

int SSsimulate2_init(SSsimulate2 *self){
	if(self==NULL)										GMERR(-11);
	memset(self,'\0',sizeof(SSsimulate2));
	return(0);
GMERRH("SSsimulate2_init",1);
}

int SSsimulate2_free(SSsimulate2 *s){
	if(s==NULL)											GMERR(-1);
	if(&GM_FreeGSLVector==NULL)							GMERR(-2);
	if(&GM_FreeGSLMatrix==NULL)							GMERR(-3);
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
						&s->eps,&s->alphaplus0,&s->at_Vm,
						&s->Zat_Vp,&s->zerosNp,
						&s->zerosNr,&s->eta_Vr,
						&s->TTat_Vm,&s->RReta_Vm,
						NULL))							GMERR(-51);
	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
						&s->y_plus_Mtp,&s->HH_C,&s->QQ_C,
						&s->alpha_plus,&s->alpha_hat_plus,
						&s->alpha_hat,&s->P_init_C,
						&s->P_init_C,NULL))				GMERR(-61);
	marray3d_free(s->K1);
	marray3d_free(s->Finv1);
	marray3d_free(s->K2);
	marray3d_free(s->Finv2);
	gsl_rng_free(s->rand);
	if(FSS2_free(s->fss1))								GMERR(-81);
	if(FSS2_free(s->fss2))								GMERR(-91);
	if(Kalman2_free(s->kalman1))						GMERR(-101);
	if(Kalman2_free(s->kalman2))						GMERR(-111);

	return(0);
GMERRH("SSsimulate2_free",1);
}

#define M_ALLOC(M,SIZE1,SIZE2) self->M = gsl_matrix_alloc(SIZE1,SIZE2)
#define V_ALLOC(V,SIZE) self->V = gsl_vector_alloc(SIZE);
#define A_ALLOC3(A,SIZE1,SIZE2,SIZE3) self->A = marray3d_alloc(SIZE1,SIZE2,SIZE3)

int SSsimulate2_construct(SSsimulate2 *self,gsl_matrix *y,gsl_vector *a_init,gsl_matrix *P_init,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,gsl_matrix *alpha_draw){
	int Nr,Nm,Np,Nt;
	if(self==NULL)										GMERR(-1);
	self->Nr = Q->d2;
	self->Nm = a_init->size;
	self->Np = y->size2;
	self->Nt = y->size1;
	Nr = self->Nr;
	Nm = self->Nm;
	Np = self->Np;
	Nt = self->Nt;

	//Matrices
	//Input
	self->y = y;
	self->P_init = P_init;
	self->alpha_draw = alpha_draw;

	//Temporary
	M_ALLOC(y_plus_Mtp,Nt,Np);
	M_ALLOC(alpha_plus,Nt,Nm);
	M_ALLOC(alpha_hat_plus,Nt,Nm);
	M_ALLOC(alpha_hat,Nt,Nm);
	M_ALLOC(v,Nt,Np);
	M_ALLOC(P_init_C,Nm,Nm);
	M_ALLOC(HH_C,Np,Np);
	M_ALLOC(QQ_C,Nr,Nr);


	//vectors
	//Input
	self->a_init = a_init;

	//Temporary
	V_ALLOC(at_Vm,Nm);
	V_ALLOC(eps,Np);
	V_ALLOC(Zat_Vp,Np);
	V_ALLOC(eta_Vr,Nr);
	V_ALLOC(TTat_Vm,Nm);
	V_ALLOC(RReta_Vm,Nm);
	V_ALLOC(alphaplus0,Nm);
	V_ALLOC(zerosNp,Np);
	V_ALLOC(zerosNr,Nr);
	gsl_vector_set_zero(self->zerosNp);
	gsl_vector_set_zero(self->zerosNr);


	//marray
	//Input
	self->H = H;
	self->T = T;
	self->R = R;
	self->Q = Q;
	self->Z = Z;

	//Temporary
	A_ALLOC3(K1,Nt,Nm,Np);
	A_ALLOC3(Finv1,Nt,Np,Np);
	A_ALLOC3(K2,Nt,Nm,Np);
	A_ALLOC3(Finv2,Nt,Np,Np);

	//Random number generator
	self->rand = gsl_rng_alloc(gsl_rng_ranlxs0);

	//Objects for filtering and smoothing
	self->kalman1 = Kalman2_New();
	self->kalman2 = Kalman2_New();
	self->fss1 = FSS2_New();
	self->fss2 = FSS2_New();
	//self->kfs1 = KFS_New();
//	self->kfs2 = KFS_New();
	printf(">>>1\n");
//	if(self->kfs1 == NULL)				printf("kfs1\n");
	if(self->y == NULL)					printf("y\n");
	if(self->a_init == NULL)			printf("a_init\n");
	if(self->P_init == NULL)			printf("P_init\n");
	if(self->Z == NULL)					printf("Z\n");
	if(self->H == NULL)					printf("H\n");
	if(self->T == NULL)					printf("T\n");
	if(self->R == NULL)					printf("R\n");
	if(self->Q == NULL)					printf("Q\n");
	if(self->v == NULL)					printf("v\n");
	if(self->K1 == NULL)				printf("K1\n");
	if(self->Finv1 == NULL)				printf("Finv1\n");
	if(self->alpha_hat == NULL)			printf("alpha_hat\n");
	printf("----------\n");

//	if(KFS_construct(self->kfs1,y,a_init,P_init,
//					Z,H,T,R,Q,self->v,self->K1,self->Finv1,self->alpha_hat))		GMERR(-11);
//	printf(">>>2\n");
//	if(KFS_construct(self->kfs2,self->y_plus_Mtp,a_init,P_init,
//					Z,H,T,R,Q,self->v,self->K2,self->Finv2,self->alpha_hat_plus))	GMERR(-11);
//	printf(">>>3\n");

	if(Kalman2_construct(self->kalman1,y,a_init,P_init,Z,H,T,R,Q,self->v,self->K1,self->Finv1))					GMERR(-11);
	if(Kalman2_construct(self->kalman2,self->y_plus_Mtp,a_init,P_init,Z,H,T,R,Q,self->v,self->K2,self->Finv2))	GMERR(-21);
	if(FSS2_construct(self->fss1,self->v,self->K1,self->Finv1,a_init,P_init,Z,T,R,Q,self->alpha_hat))			GMERR(-31);
	if(FSS2_construct(self->fss2,self->v,self->K2,self->Finv2,a_init,P_init,Z,T,R,Q,self->alpha_hat_plus))		GMERR(-31);

	return(0);
GMERRH("SSsimulate2_construct",1);
}

#undef M_ALLOC
#undef V_ALLOC
#undef A_ALLOC3

int SSsimulate2_Simulate(gsl_matrix *y,gsl_vector *a_init,gsl_matrix *P_init,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,gsl_matrix *alpha_draw){
	SSsimulate2 *self =NULL;
	if((self=SSsimulate2_New())==NULL)						GMERR(-11);
	if(SSsimulate2_construct(self,y,a_init,P_init,Z,H,T,R,Q,alpha_draw))	GMERR(-21);
	if(_SSsimulate2_checkSizes(self))						GMERR(-31);
	if(SSsimulate2_operations(self))							GMERR(-41);
	if(SSsimulate2_free(self))								GMERR(-51);
	return(0);
GMERRH("SSsimulate2_Simulate",1);
}

/*Checks the sizes of the inputs and outputs to make sure
 *that they match relative to y and only used in the construct
 *method
*/

static int _SSsimulate2_checkSizes(SSsimulate2 *self){
	int Nt,Nr,Nm,Np;
	if(self==NULL)									GMERR(-1);
	Np=self->Np;Nr=self->Nr;Nm=self->Nm;Nt=self->Nt;

	if(marray3d_checkYZ(self->H,Np,Np))				GMERR(-11);		
	if(marray3d_checkYZ(self->T,Nm,Nm))				GMERR(-21);
	if(marray3d_checkYZ(self->Q,Nr,Nr))				GMERR(-31);
	if(marray3d_checkYZ(self->R,Nm,Nr))				GMERR(-41);
	if(marray3d_checkYZ(self->Z,Np,Nm))				GMERR(-51);
	if(matrixCheckSize(self->P_init,Nm,Nm))			GMERR(-61);
	if(matrixCheckSize(self->y,Nt,Np))				GMERR(-71);
	if(marray3d_checkX1orT(self->H,Nt))				GMERR(-81);
	if(marray3d_checkX1orT(self->R,Nt))				GMERR(-91);
	if(marray3d_checkX1orT(self->T,Nt))				GMERR(-101);
	if(marray3d_checkX1orT(self->Q,Nt))				GMERR(-111);
	if(marray3d_checkX1orT(self->Z,Nt))				GMERR(-121);

	return(0);
GMERRH("_SSsimulate2_checkSizes",1);
}


/*Zeros out any of the outputs to make sure that the stuff comes from
 *this run
*/
static int _SSsimulate2_zeroOutputs(SSsimulate2 *self){
	if(self==NULL)											GMERR(-1);
	if(self->alpha_draw==NULL)								GMERR(-2);
	gsl_matrix_set_zero(self->alpha_draw);
	return(0);
GMERRH("_SSsimulate2_zeroOutputs",1);
}

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

int SSsimulate2_operations(SSsimulate2 *self){
	int t;
	int Tt,Rt,Ht,Qt,Zt;
	int Nt,Np,Nr,Nm;
	gsl_matrix *TT,*HH,*RR,*QQ,*ZZ;
	if(self==NULL)															GMERR(-1);
	Nt = self->Nt;
	Np = self->Np;
	Nr = self->Nr;
	Nm = self->Nm;
	//gsl_vector *alphaEnd = gsl_vector_alloc(Np);
	//Used to determine if we must update the transition matrices
	Ht = ((self->H->d1)>1);
	Tt = ((self->T->d1)>1);
	Rt = ((self->R->d1)>1);
	Qt = ((self->Q->d1)>1);
	Zt = ((self->Z->d1)>1);

	if(_SSsimulate2_zeroOutputs(self))										GMERR(-11);

	//alpha_plus[0] = np.random.multivariate_normal(a_init,P_init)
	//if(rmvnorm(self->rand,Nm,self->a_init,self->P_init,self->alphaplus0))	GMERR(-21);
	gsl_matrix_memcpy(self->P_init_C,self->P_init);
	gsl_linalg_cholesky_decomp1(self->P_init_C);
	if(gsl_ran_multivariate_gaussian(self->rand,self->a_init,self->P_init_C,self->alphaplus0))	GMERR(-21);

	//gsl_matrix_get_row(self->alphaplus0,self->alpha0Samps,self->a);
	if(M_SET_ROW(alpha_plus,alphaplus0,0))									GMERR(-31);
	
	//set initial values for TT,RT,HT,QT
	A_ROW(T,TT,0);
	A_ROW(H,HH,0);
	A_ROW(R,RR,0);
	A_ROW(Q,QQ,0);
	A_ROW(Z,ZZ,0);

	for(t=0;t<Nt;t++){
		//Updating transition matrices if necessary
		if(Tt) A_ROW(T,TT,t);
		if(Rt) A_ROW(R,RR,t);
		if(Ht) A_ROW(H,HH,t);
		if(Qt) A_ROW(Q,QQ,t);
		if(Zt) A_ROW(Z,ZZ,t);

		//eps = np.random.multivariate_normal(np.zeros(Np),HH)
		//if(rmvnorm(self->rand,Np,self->zerosNp,HH,self->eps))				GMERR(-41);
		gsl_matrix_memcpy(self->HH_C,HH);
		gsl_linalg_cholesky_decomp1(self->HH_C);
		if(gsl_ran_multivariate_gaussian(self->rand,self->zerosNp,self->HH_C,self->eps))	GMERR(-41);

		//y_plus[t] = ZZ.dot(alpha_plus[t]) + eps
		if(M_GET_ROW(alpha_plus,at_Vm,t))									GMERR(-51);
		gsl_blas_dgemv(CblasNoTrans,1.0,ZZ,self->at_Vm,0.0,self->Zat_Vp);
		//if(DDOT_MV(ZZ,at_Vm,Zat_Vp))										GMERR(-61);
		if(V_ADD(Zat_Vp,eps))												GMERR(-71);
		if(M_SET_ROW(y_plus_Mtp,Zat_Vp,t))									GMERR(-81);

		if(t+1<Nt){
			//eta = np.random.multivariate_normal(np.zeros(Nr),QQ)
			//if(rmvnorm(self->rand,Nr,self->zerosNr,QQ,self->eta_Vr))		GMERR(-171);
			gsl_matrix_memcpy(self->QQ_C,QQ);
			gsl_linalg_cholesky_decomp1(self->QQ_C);
			if(gsl_ran_multivariate_gaussian(self->rand,self->zerosNr,self->QQ_C,self->eta_Vr))	GMERR(-171);

			//alpha_plus[t+1] = TT.dot(alpha_plus[t]) + RR.dot(eta)
			//if(DDOT_MV(TT,at_Vm,TTat_Vm))									GMERR(-181);
			gsl_blas_dgemv(CblasNoTrans,1.0,TT,self->at_Vm,0.0,self->TTat_Vm);
			//if(DDOT_MV(RR,eta_Vr,RReta_Vm))									GMERR(-191);
			gsl_blas_dgemv(CblasNoTrans,1.0,RR,self->eta_Vr,0.0,self->RReta_Vm);
			if(V_ADD(TTat_Vm,RReta_Vm))										GMERR(-201);
			if(M_SET_ROW(alpha_plus,TTat_Vm,(t+1)))							GMERR(-211);
		}
	}

	if(self->interleaved){
		printf("No methods for interleaved\n");
		GMERR(-221);
	}else{
	//	if(KFS_operations(self->kfs1))				GMERR(-221);
	//	if(KFS_operations(self->kfs2))				GMERR(-231);
		if(SSsimulate2_runParallel(self))			GMERR(-221);

		if(Kalman2_operations(self->kalman1))								GMERR(-221);
		if(FSS2_operations(self->fss1))										GMERR(-231);
		if(Kalman2_operations(self->kalman2))								GMERR(-241);
		if(FSS2_operations(self->fss2))										GMERR(-251);

		if(M_ADD(alpha_draw,alpha_plus))									GMERR(-261);
		if(M_SUB(alpha_draw,alpha_hat_plus))								GMERR(-271);

		if(M_ADD(alpha_draw,alpha_hat))										GMERR(-281);

	}
	return(0);
GMERRH("SSsimulate2_operations",1);
}

static void *_run1(void *_self){
	SSsimulate2 *self = (SSsimulate2 *)_self;
//	if(KFS_operations(self->kfs1))	fprintf(stderr,"Error pthread 1\n");
	return(NULL);
}

static void *_run2(void *_self){
	SSsimulate2 *self = (SSsimulate2 *)_self;
//	if(KFS_operations(self->kfs2))	fprintf(stderr,"Error pthread 2\n");	
	return(NULL);
}

int SSsimulate2_runParallel(SSsimulate2 *self){
	pthread_t myThread1,myThread2;
	if(pthread_create(&myThread1,NULL,&_run1,(void *)self))		GMERR(-11);
	if(pthread_create(&myThread2,NULL,&_run2,(void *)self))		GMERR(-21);
	pthread_join(myThread1,NULL);
	pthread_join(myThread2,NULL);
	return(0);
GMERRH("SSsimulate2_runParallel",1);
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

int SSsimulate2_writeAlphaDraw(SSsimulate2 *self,char *name){
	FILE *fp;
	if(self==NULL)								GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)								GMERR(-3);
	if(gsl_matrix_fprintf(fp,self->alpha_draw,"%.15f"))	GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("SSsimulate2_writeAlphaDraw",1);
}

int SSsimulate2_writeOutputs(SSsimulate2 *self,char *baseName){
	FILE *fp;
	char strA[16] = "_AlphaDraw.txt";
	char outA[1024];
	sprintf(outA,"%s%s",baseName,strA);
	if(self==NULL)								GMERR(-1);
	if(SSsimulate2_writeAlphaDraw(self,outA))	GMERR(-11);
	return(0);
GMERRH("SSsimulate2_writeOutputs",1);
}

