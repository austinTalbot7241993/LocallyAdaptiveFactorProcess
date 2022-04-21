#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_exp.h>
#include "gmlib.h"
#include "matrix_utils.h"
#include "NGPmcmc.h"
#include "mdarray.h"
#include "KalmanFilter2.h"
#include "fastStateSmoother2.h"
#include "SSsimulate2.h"
#include "NGPtools2.h"

/*Code by Austin Bryan Talbot
 *		  a.talbot@verizon.net
 *
 *Description:
 *
 *Completed:
 *
 *Revision History:
 *
*/

#define NGPmcmc_MIN(A,B) ({\
	typeof(A) _A_TEMP_;\
	typeof(B) _B_TEMP_;\
	_A_TEMP_ = (A);\
	_B_TEMP_ = (B);\
	_A_TEMP_ = _A_TEMP_ > _B_TEMP_ ? _B_TEMP_ : _A_TEMP_;\
	})



//This is used for error handling by macros in gmlib
static int _NGPmcmc_checkSizes(gsl_matrix *y,gsl_vector *tobs,int Niter,gsl_vector *sigU,gsl_vector *sigA,double sigEps,double sigMu,double sigAlph,double a,double b,marray3d *th,gsl_matrix *sig);


static int _NGPmcmc_zeroOutputs(NGPmcmc *self);

//Sampling methods
int NGPmcmc_drawState(NGPmcmc *self);
int NGPmcmc_drawSigEps(NGPmcmc *self);
int NGPmcmc_calc_acceptance_prob(NGPmcmc *self);
int NGPmcmc_drawSigs(NGPmcmc *self);

/******************************************************
* This section defines the methods for instantiating, *
* creating, and freeing the object 					  *
******************************************************/


NGPmcmc * NGPmcmc_New(){
	NGPmcmc *self;
	self = (NGPmcmc *)GM_Malloc(sizeof(NGPmcmc));
	if(self==NULL)							GMERR(-11);
	if(NGPmcmc_init(self))					GMERR(-21);
	return(self);
GMERRH("NGPmcmc_New",NULL);
}

int NGPmcmc_init(NGPmcmc *self){
	if(self==NULL)							GMERR(-1);
	memset(self,'\0',sizeof(NGPmcmc));

	return(0);
GMERRH("NGPmcmc_init",1);
}

int NGPmcmc_free(NGPmcmc *s){
	if(s==NULL)								GMERR(-1);
	if(&GM_FreeGSLVector==NULL)				GMERR(-2);
	if(&GM_FreeGSLMatrix==NULL)				GMERR(-3);

	gsl_vector_free(s->y_vec);
	gsl_vector_free(s->delta);

	/*************
	* draw_state *
	*************/
	marray3d_free(s->Z_s);
	marray3d_free(s->H_s);
	marray3d_free(s->T_s);
	marray3d_free(s->R_s);
	marray3d_free(s->Q_s);
	gsl_vector_free(s->a0_s);
	gsl_matrix_free(s->P0_s);
	if(SSsimulate2_free(s->drawState))				GMERR(-11);

	/**************
	* draw_sigEps *
	**************/
	gsl_vector_free(s->U);

	/***********************
	* calc_acceptance_prob *
	***********************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->v_t,&s->vstar_t,&s->v_03,
				&s->thetastar_t,&s->theta_t,
				&s->v_02,&s->work3,&s->work2,
				&s->Gthetastar,NULL))	GMERR(-21);

	gsl_matrix_free(s->HtildeT);

	/************
	* draw_sigs *
	************/
	gsl_vector_free(s->sigUstar);
	gsl_vector_free(s->sigAstar);

	gsl_matrix_free(s->theta);
	gsl_matrix_free(s->thetastar);

	gsl_vector_free(s->Ustar);
	gsl_vector_free(s->DUstar);
	gsl_vector_free(s->Astar);

	//Regular
	marray3d_free(s->Zdummy);
	marray3d_free(s->Hdummy);
	marray3d_free(s->G);
	marray3d_free(s->H);
	marray3d_free(s->W);
	gsl_vector_free(s->a0dummy);
	gsl_matrix_free(s->P0dummy);

	//tilde
	marray3d_free(s->Ztdummy);
	marray3d_free(s->Htdummy);
	marray3d_free(s->Gtilde);
	marray3d_free(s->Htilde);
	marray3d_free(s->Wtilde);
	gsl_vector_free(s->a0tildedummy);
	gsl_matrix_free(s->P0tildedummy);

	//Star
	marray3d_free(s->Zstardummy);
	marray3d_free(s->Hstardummy);
	marray3d_free(s->Gstardummy);
	marray3d_free(s->Tstardummy);
	marray3d_free(s->Wstar);
	gsl_vector_free(s->a0stardummy);
	gsl_matrix_free(s->P0stardummy);

	//tilde star
	marray3d_free(s->Ztildestardummy);
	marray3d_free(s->Htildestardummy);
	marray3d_free(s->Gtildestardummy);
	marray3d_free(s->Ttildestardummy);
	marray3d_free(s->Wtildestar);
	gsl_vector_free(s->a0tildestardummy);
	gsl_matrix_free(s->P0tildestardummy);

	marray3d_free(s->Z_star);
	marray3d_free(s->H_star);
	marray3d_free(s->T_star);
	marray3d_free(s->R_star);
	marray3d_free(s->Q_star);
	gsl_vector_free(s->a0_star);
	gsl_matrix_free(s->P0_star);
	if(SSsimulate2_free(s->drawThetaStar))	GMERR(-41);

	gsl_rng_free(s->rand);

	return(0);
GMERRH("NGPmcmc_free",1);
}

/*******************************************
* Any specific auxillary memory allocation *
* methods are here 						   *
*******************************************/

/***************************************
* Auxilary method to check the input   *
* sizes before constructing the object *
***************************************/
static int _NGPmcmc_checkSizes(gsl_matrix *y,gsl_vector *tobs,int Niter,gsl_vector *sigU,gsl_vector *sigA,double sigEps,double sigMu,double sigAlph,double a,double b,marray3d *th,gsl_matrix *sig){
	int Nt;
	if((y==NULL)||(tobs==NULL)||(sigU==NULL)||
		(sigA==NULL)||(th==NULL)||(sig==NULL))	GMERR(-1);
	Nt = y->size1;
	
	//Size checking
	if(y->size2!=1)								GMERR(-11);
	if(tobs->size!=Nt)							GMERR(-21);
	if(sigU->size!=1)							GMERR(-31);
	if(sigA->size!=1)							GMERR(-41);
	if(marray3d_checkSizes(th,Niter,Nt,3))		GMERR(-51);
	if(matrixCheckSize(sig,Niter,3))			GMERR(-61);

	//Reasonable answers
	if(sigEps<=0)								GMERR(-71);
	if(sigAlph<=0)								GMERR(-81);
	if(sigMu<=0)								GMERR(-91);
	if(b<=0)									GMERR(-101);
	if(a<=0)									GMERR(-111);
	if(positiveVector(sigU))					GMERR(-121);
	if(positiveVector(sigA))					GMERR(-131);

	return(0);
GMERRH("_NGPmcmc_checkSizes",1);
}

/*****************************************
* Zeros the outputs. I use this whenever *
* I change a parameter so that if I 	 *
* forget to rerun it I dont use wrong 	 *
* results 								 *
*****************************************/

static int _NGPmcmc_zeroOutputs(NGPmcmc *self){
	if(self==NULL)						GMERR(-1);
	if(marray3d_set_zero(self->th))		GMERR(-11);
	gsl_matrix_set_zero(self->sig);

	return(0);
GMERRH("_NGPmcmc_zeroOutputs",1);
}

//Memory allocation macros
#define V_ALLOC(V,SIZE) self->V = gsl_vector_alloc(SIZE)
#define M_ALLOC(M,SIZE1,SIZE2) self->M = gsl_matrix_alloc(SIZE1,SIZE2)
#define A_ALLOC3(A,SIZE1,SIZE2,SIZE3) self->A = marray3d_alloc(SIZE1,SIZE2,SIZE3)

int NGPmcmc_construct(NGPmcmc *self,gsl_matrix *y,gsl_vector *tobs,int Niter,gsl_vector *sigU,gsl_vector *sigA,double sigEps,double sigMu,double sigAlph,double a,double b,marray3d *th,gsl_matrix *sig){
	int Nt,t;
	double temp;
	if(_NGPmcmc_checkSizes(y,tobs,Niter,sigU,sigA,sigEps,sigMu,sigAlph,a,b,th,sig))	GMERR(-11);

	/******************
	*******************
	** method inputs **
	*******************
	******************/
	self->y = y;
	self->tobs = tobs;
	self->Niter = Niter;
	self->sigU = sigU;
	self->sigA = sigA;
	self->sigEps = sigEps;
	self->sigMu = sigMu;
	self->sigAlph = sigAlph;
	self->a = a;
	self->Nt = y->size1;
	Nt = self->Nt;
	self->burn = 1000;

	/*******************
	********************
	** method outputs **
	********************
	*******************/
	self->th = th;
	self->sig = sig;

	/*********
	* priors *
	*********/
	self->a = a;
	self->b = b;
	self->sigMu = sigMu;
	self->sigAlph = sigAlph;

	/*******************************
	********************************
	** method temporary variables **
	********************************
	*******************************/
	V_ALLOC(y_vec,Nt);
	V_ALLOC(delta,Nt);
	for(t=0;t<(Nt-1);t++){
		temp = gsl_vector_get(tobs,t+1)-gsl_vector_get(tobs,t);
		gsl_vector_set(self->delta,t,temp);
	}
	gsl_vector_set(self->delta,Nt-1,gsl_vector_get(self->delta,Nt-2));

	gsl_matrix_get_col(self->y_vec,y,0);

	/************
	* draw_sigs *
	************/

	V_ALLOC(sigUstar,1);
	V_ALLOC(sigAstar,1);

	M_ALLOC(theta,Nt,3);
	M_ALLOC(thetastar,Nt,3);

	V_ALLOC(Ustar,Nt);
	V_ALLOC(DUstar,Nt);
	V_ALLOC(Astar,Nt);

	//Regular
	A_ALLOC3(Zdummy,1,1,3);
	A_ALLOC3(Hdummy,1,1,1);
	A_ALLOC3(G,Nt,3,3);
	A_ALLOC3(H,1,3,3);
	A_ALLOC3(W,Nt,3,3);
	V_ALLOC(a0dummy,3);
	M_ALLOC(P0dummy,3,3);

	//Tilde
	A_ALLOC3(Ztdummy,1,1,3);
	A_ALLOC3(Htdummy,1,1,1);
	A_ALLOC3(Gtilde,Nt,3,3);
	A_ALLOC3(Htilde,1,3,2);
	A_ALLOC3(Wtilde,Nt,2,2);
	V_ALLOC(a0tildedummy,3);
	M_ALLOC(P0tildedummy,3,3);

	//Star
	A_ALLOC3(Zstardummy,1,1,3);
	A_ALLOC3(Hstardummy,1,1,1);
	A_ALLOC3(Gstardummy,Nt,3,3);
	A_ALLOC3(Tstardummy,1,3,3);
	A_ALLOC3(Wstar,Nt,3,3);
	V_ALLOC(a0stardummy,3);
	M_ALLOC(P0stardummy,3,3);

	//Tilde star
	A_ALLOC3(Ztildestardummy,1,1,3);
	A_ALLOC3(Htildestardummy,1,1,1);
	A_ALLOC3(Gtildestardummy,Nt,3,3);
	A_ALLOC3(Ttildestardummy,1,3,2);
	A_ALLOC3(Wtildestar,Nt,2,2);
	V_ALLOC(a0tildestardummy,3);
	M_ALLOC(P0tildestardummy,3,3);

	A_ALLOC3(Z_star,1,1,3);
	A_ALLOC3(H_star,1,1,1);
	A_ALLOC3(T_star,Nt,3,3);
	A_ALLOC3(R_star,1,3,2);
	A_ALLOC3(Q_star,Nt,2,2);
	V_ALLOC(a0_star,3);
	M_ALLOC(P0_star,3,3);

	/**************
	* draw_sigEps *
	**************/
	V_ALLOC(U,Nt);

	/***********************
	* calc_acceptance_prob *
	***********************/
	self->prob = 0;
	V_ALLOC(v_t,3);
	V_ALLOC(vstar_t,2);
	V_ALLOC(v_03,3);
	gsl_vector_set_zero(self->v_03);
	V_ALLOC(v_02,2);
	gsl_vector_set_zero(self->v_02);
	V_ALLOC(work3,3);
	V_ALLOC(work2,2);
	V_ALLOC(theta_t,3);
	V_ALLOC(thetastar_t,3);
	V_ALLOC(Gthetastar,3);

	M_ALLOC(HtildeT,2,3);

	/*************
	* draw_state *
	*************/

	self->drawState = SSsimulate2_New();
	self->drawThetaStar = SSsimulate2_New();

	V_ALLOC(a0_s,3);
	M_ALLOC(P0_s,3,3);
	A_ALLOC3(Z_s,1,1,3);
	A_ALLOC3(H_s,1,1,1);
	A_ALLOC3(T_s,Nt,3,3);
	A_ALLOC3(R_s,1,3,3);
	A_ALLOC3(Q_s,Nt,3,3);

	if(SSsimulate2_construct(self->drawState,y,self->a0_s,self->P0_s,
							self->Z_s,self->H_s,self->T_s,self->R_s,
							self->Q_s,self->theta))		GMERR(-11);

	if(SSsimulate2_construct(self->drawThetaStar,y,self->a0_star,
							self->P0_star,self->Z_star,self->H_star,
							self->T_star,self->R_star,self->Q_star,
							self->thetastar))			GMERR(-21);


	/**************************
	* random number generator *
	**************************/
	self->rand = gsl_rng_alloc(gsl_rng_mt19937);

	return(0);
GMERRH("NGPmcmc_construct",1);
}

#undef V_ALLOC
#undef M_ALLOC
#undef A_ALLOC3

/*********************************************************
* Here is the method as you would see it in conventional *
* programming styles.									 *
*********************************************************/

//Mathematical macros
#define A_GET_ROW(A,M,T) marray3d_get_X(self->A,self->M,T)
#define A_SET_ROW(A,M,T) marray3d_set_X(self->A,self->M,T)
#define DDOT_MV(M,V,R) gsl_blas_dgemv(CblasNoTrans,1.0,self->M,self->V,0.0,self->R)
#define DDOT_VV(UU,VV,RR) gsl_blas_ddot(self->UU,self->VV,RR)
#define DDOT_MM(MM,NN,RR) gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,self->MM,self->NN,0.0,self->RR)
#define M_GET_ROW(M,V,T) gsl_matrix_get_row(self->V,self->M,T)
#define M_SET_ROW(M,V,T) gsl_matrix_set_row(self->M,T,self->V)
#define M_ADD(MM,NN) gsl_matrix_add(self->MM,self->NN)
#define V_SUB(VV,UU) gsl_vector_sub(self->VV,self->UU)
#define V_ADD(VV,UU) gsl_vector_add(self->VV,self->UU)
#define M_SUB(MM,NN) gsl_matrix_sub(self->MM,self->NN)
#define DDOT_VVOUTER(VV,MM) vector_outer(self->VV,self->MM)
#define A_ROW(A,M,T) M = self->A->matrixList[T]

int NGPmcmc_operations(NGPmcmc *self){
	int t,burn,iter;
	if(self==NULL)					GMERR(-1);
	burn = self->burn;
	iter = self->Niter;

	//th = Array(Float64,Nm,Nt,Niter+burnin)
	//sig = Array(Float64,3,Niter+burnin)
	if(_NGPmcmc_zeroOutputs(self))	GMERR(-11);

	for(t=0;t<burn;t++){
		if((t%1000==0))	printf("Burnin %d\n",t);
		//theta = draw_state(sigEps,sigU,sigA)
		if(NGPmcmc_drawState(self))	GMERR(-111);

		//sigEps = draw_sigEps(theta[1,:])
		if(NGPmcmc_drawSigEps(self))GMERR(-121);

		//sigU,sigA = draw_sigs(theta,sigEps,sigU,sigA)
		if(NGPmcmc_drawSigs(self))	GMERR(-131);

	}
	for(t=0;t<iter;t++){
		if((t%1000==0))	printf("Iteration %d\n",t);
		//theta = draw_state(sigEps,sigU,sigA)
		if(NGPmcmc_drawState(self))	GMERR(-111);

		//sigEps = draw_sigEps(theta[1,:])
		if(NGPmcmc_drawSigEps(self))GMERR(-121);

		//sigU,sigA = draw_sigs(theta,sigEps,sigU,sigA)
		if(NGPmcmc_drawSigs(self))	GMERR(-131);

		//th[:,:,idx] = theta
		if(marray3d_set_X(self->th,self->theta,t))	GMERR(-241);

		//sig[:,idx] = [sigEps,sigU,sigA]
		gsl_matrix_set(self->sig,t,0,self->sigEps);
		gsl_matrix_set(self->sig,t,1,gsl_vector_get(self->sigU,0));
		gsl_matrix_set(self->sig,t,2,gsl_vector_get(self->sigA,0));
	}
	
	return(0);
GMERRH("NGPmcmc_operations",1);
}

int NGPmcmc_drawState(NGPmcmc *self){
	if(self==NULL)				GMERR(-1);
	//Z,H,T,R,Q,a-init,Pinit = assemble_matrices(dims,delta,sigeps,sigU,sigA,sigmu,sigalpha,approx)
	if(NGP_assemble_matrices(1,1,self->delta,&self->sigEps,
			self->sigU,self->sigA,&self->sigMu,&self->sigAlph,
			0,self->Z_s,self->H_s,self->T_s,self->R_s,
			self->Q_s,self->a0_s,self->P0_s))	GMERR(-11);

	//NGPtools.sample(y,1,Np,length(tobs),delta,sigEps,sigU,sigA,sigmu,sigalph,approx=approx)[:,:,1]
	if(SSsimulate2_operations(self->drawState))	GMERR(-21);

	return(0);
GMERRH("NGPmcmc_drawState",1);
}

//checked
int NGPmcmc_drawSigEps(NGPmcmc *self){
	double a_eff=0,b_eff=0,temp=0;
	if(self==NULL)				GMERR(-1);
	
	gsl_matrix_get_col(self->U,self->theta,0);

	//a_eff = a + (Nt/2)
	a_eff = self->a + 0.5*(double)self->Nt;

	//b_eff = b + 0.5 * sum((y-U).^2)
	if(vector_diff2(self->y_vec,self->U,&b_eff))		GMERR(-11);
	b_eff = 0.5*b_eff;
	b_eff = b_eff + self->b;

	//return(sqrt(rand(InverseGamma(a_eff,b_eff))))
	temp = sqrt(1.0/gsl_ran_gamma(self->rand,a_eff,1.0/b_eff));
	self->sigEps = temp;

	return(0);
GMERRH("NGPmcmc_drawSigEps",1);
}

//checked
int NGPmcmc_calc_acceptance_prob(NGPmcmc *self){
	int Nt,updateG,updateW,updateWstar,updateWtilde,updateWtildestar;
	gsl_matrix *myGt,*myHtilde,*myWs,*myWw,*myWt,*myWts;
	double pi,temp;
	int t;
	if(self==NULL)				GMERR(-1);

	//Variables that need to be updated
	updateG = (self->G->d1>1);
	updateW = (self->W->d1>1);
	updateWstar = (self->Wstar->d1>1);
	updateWtilde = (self->Wtilde->d1>1);
	updateWtildestar = (self->Wtildestar->d1>1);

	Nt = self->Nt;

	//Note that computing the cholesky decomposition destroys the original matrix
	//I save time from not copying memory and note that we never use these after
	//this method is run in draw sigs

	gsl_linalg_cholesky_decomp1(self->W->matrixList[0]);
	gsl_linalg_cholesky_decomp1(self->Wstar->matrixList[0]);
	gsl_linalg_cholesky_decomp1(self->Wtilde->matrixList[0]);
	gsl_linalg_cholesky_decomp1(self->Wtildestar->matrixList[0]);

	for(t=1;t<Nt;t++){
		if(updateW)				gsl_linalg_cholesky_decomp1(self->W->matrixList[t]);
		if(updateWstar)			gsl_linalg_cholesky_decomp1(self->Wstar->matrixList[t]);
		if(updateWtilde)		gsl_linalg_cholesky_decomp1(self->Wtilde->matrixList[t]);
		if(updateWtildestar)	gsl_linalg_cholesky_decomp1(self->Wtildestar->matrixList[t]);
	}

	//Initializing all the temporary pointers
	A_ROW(G,myGt,0);
	A_ROW(Htilde,myHtilde,0);
	A_ROW(Wstar,myWs,0);
	A_ROW(W,myWw,0);
	A_ROW(Wtilde,myWt,0);
	A_ROW(Wtildestar,myWts,0);

	for(t=0;t<(Nt-1);t++){
		//local G_t = ndims(G) < 3 ? G : slice(G,:,:,t)
		if(updateG)		A_ROW(G,myGt,t);

		//v[:,t] = theta[:,t+1] - G_t *Theta[:,t]
		if(M_GET_ROW(theta,theta_t,t))		GMERR(-11);
		gsl_blas_dgemv(CblasNoTrans,1.0,myGt,self->theta_t,0.0,self->v_t);
		gsl_vector_scale(self->v_t,-1.0);
		if(M_GET_ROW(theta,theta_t,t+1))	GMERR(-21);
		gsl_vector_add(self->v_t,self->theta_t);

		//v_star[:,t] = Htilde' * (Thetastar[:,t+1]-G_t * thetastar[:,t])
		if(M_GET_ROW(thetastar,thetastar_t,t))	GMERR(-31);
		gsl_blas_dgemv(CblasNoTrans,1.0,myGt,self->thetastar_t,0.0,self->Gthetastar);
		gsl_vector_scale(self->Gthetastar,-1.0);
		if(M_GET_ROW(thetastar,thetastar_t,t+1))GMERR(-41);
		gsl_vector_add(self->Gthetastar,self->thetastar_t);
		if(gTranspose(myHtilde,self->HtildeT))					GMERR(-999);
		gsl_blas_dgemv(CblasNoTrans,1.0,self->HtildeT,self->Gthetastar,0.0,self->vstar_t);
		
//		printGSLVectorT(self->vstar_t);

		if(updateW)				A_ROW(W,myWw,t);
		if(updateWstar)			A_ROW(Wstar,myWs,t);
		if(updateWtilde)		A_ROW(Wtilde,myWt,t);
		if(updateWtildestar)	A_ROW(Wtildestar,myWts,t);

		//pi += stats.multivaraite_normal.logpdf(x=v[t],mean=zeros,cov=Wstar[t])
		gsl_ran_multivariate_gaussian_log_pdf(self->v_t,self->v_03,myWs,&temp,self->work3);
		pi += temp;

		//pi -= stats.multivaraite_normal.logpdf(x=v[t],mean=None,cov=W[t])
		gsl_ran_multivariate_gaussian_log_pdf(self->v_t,self->v_03,myWw,&temp,self->work3);
		pi -= temp;

		//pi += stats.multivaraite_normal.logpdf(x=v_star[t],mean=zeros,cov=Wtilde[t])
		gsl_ran_multivariate_gaussian_log_pdf(self->vstar_t,self->v_02,myWt,&temp,self->work2);
		pi += temp;

		//pi -= stats.multivaraite_normal.logpdf(x=v_star[t],mean=zeros,cov=Wtilde_star[t])
		gsl_ran_multivariate_gaussian_log_pdf(self->vstar_t,self->v_02,myWts,&temp,self->work2);
		pi -= temp;
	}

	gsl_ran_multivariate_gaussian_log_pdf(self->v_03,self->v_03,myWs,&temp,self->work3);
	pi += temp;
	gsl_ran_multivariate_gaussian_log_pdf(self->v_03,self->v_03,myWw,&temp,self->work3);
	pi -= temp;
	gsl_ran_multivariate_gaussian_log_pdf(self->v_02,self->v_02,myWt,&temp,self->work2);
	pi += temp;
	gsl_ran_multivariate_gaussian_log_pdf(self->v_02,self->v_02,myWts,&temp,self->work2);
	pi -= temp;

	if(pi>1000){
		self->prob=1;
	}else if(pi<-30){
		self->prob=0;	
	}else{
		pi = gsl_sf_exp(pi);
		pi = NGPmcmc_MIN(pi,1);
		self->prob = pi;
	}
	return(0);
GMERRH("NGPmcmc_calc_acceptance_prob",1);
}

int NGPmcmc_drawSigs(NGPmcmc *self){
	int Nt,t;
	double R2U,R2A,temp,d;
	double a_eff,bUstar,bAstar,b_eff,SU,SA,u;
	if(self==NULL)				GMERR(-1);
	Nt = self->Nt;

	//Thetastar = draw_state(sigEps,sigU,sigA,approx=true)
	if(NGP_assemble_matrices(1,1,self->delta,&self->sigEps,self->sigU,self->sigA,&self->sigMu,&self->sigAlph,1,self->Z_star,self->H_star,self->T_star,self->R_star,self->Q_star,self->a0_star,self->P0_star))	GMERR(-11);
	if(SSsimulate2_operations(self->drawThetaStar))				GMERR(-21);

	//Ustar = Thetastar[1,:][:]
	gsl_matrix_get_col(self->Ustar,self->thetastar,0);

	//DUstar = Thetastar[2,:][:]
	gsl_matrix_get_col(self->DUstar,self->thetastar,1);

	//Astar = Thetastar[3,:][:]
	gsl_matrix_get_col(self->Astar,self->thetastar,2);

	//d = delta *ones(Nt-1)
	//R2U = sum((diff(DUstar)-Astar[1:(end-1)].*d).^2./d)
	//R2A = sum(diff(Astar).^2./d)
	for(t=0;t<Nt-1;t++){
		d = gsl_vector_get(self->delta,t);
		temp = (gsl_vector_get(self->DUstar,t+1)-gsl_vector_get(self->DUstar,t)-gsl_vector_get(self->Astar,t)*d);
		temp = temp*temp/d;
		R2U += temp;

		temp = (gsl_vector_get(self->Astar,t+1)-gsl_vector_get(self->Astar,t));
		temp = temp*temp/d;
		R2A += temp;
	}

	//a_eff = a + (Nt/2)
	a_eff = self->a + 0.5*(double)Nt;
	bUstar = self->b + 0.5*R2U;
	bAstar = self->b + 0.5*R2A;

	
	//sigUstar = sqrt(rand(InverseGamma(a_eff,b+R2U/2)))
	SU = sqrt(1.0/gsl_ran_gamma(self->rand,a_eff,1.0/bUstar));
	gsl_vector_set(self->sigUstar,0,SU);

	//sigAstar = sqrt(rand(InverseGamma(a_eff,b+R2A/2)))
	SA = sqrt(1.0/gsl_ran_gamma(self->rand,a_eff,1.0/bAstar));
	gsl_vector_set(self->sigAstar,0,SA);

	//_,_,G,H,W,_,_ = NGPtools.assemble_matrices(Np,delta,sigEps,sigU,sigA,sigmu,sigalph)
	if(NGP_assemble_matrices(1,1,self->delta,&self->sigEps,self->sigU,self->sigA,
	&self->sigMu,&self->sigAlph,0,self->Zdummy,self->Hdummy,self->G,self->H,
	self->W,self->a0dummy,self->P0dummy))	GMERR(-31);

	//_,_,Gtilde,Htilde,Wtilde,_,_ = NGPtools.assemble_matrices(Np,delta,sigEps,sigU,sigA,sigmu,sigalph,approx=true)

	if(NGP_assemble_matrices(1,1,self->delta,&self->sigEps,self->sigU,self->sigA,
	&self->sigMu,&self->sigAlph,1,self->Ztdummy,self->Htdummy,self->Gtilde,
	self->Htilde,self->Wtilde,self->a0tildedummy,self->P0tildedummy))	GMERR(-41);

	//_,_,_,_,Wstar,_,_ = NGPtools.assemble_matrices(Np,delta,sigEps,sigUstar,sigAstar,sigmu,sigalph)

	if(NGP_assemble_matrices(1,1,self->delta,&self->sigEps,self->sigUstar,
	self->sigAstar,&self->sigMu,&self->sigAlph,0,self->Zstardummy,self->Hstardummy,
	self->Gstardummy,self->Tstardummy,self->Wstar,self->a0stardummy,self->P0stardummy))	GMERR(-51);

	//_,_,_,_,Wtilde_star,_,_ = NGPtools.assemble_matrices(Np,delta,sigEps,sigUstar,sigAstar,sigmu,sigalph,approx=true)
	if(NGP_assemble_matrices(1,1,self->delta,&self->sigEps,self->sigU,self->sigA,&self->sigMu,
	&self->sigAlph,1,self->Ztildestardummy,self->Htildestardummy,self->Gtildestardummy,self->Ttildestardummy,self->Wtildestar,self->a0tildestardummy,self->P0tildestardummy))	GMERR(-11);

	//acc_prob = calc_acceptance_prob(theta,thetastar,G,Gtilde,W,Wstar,Wtilde,Wtilde_star,Htilde)
	if(NGPmcmc_calc_acceptance_prob(self))					GMERR(-999);

	//if rand() < acc_prob
		//return sigUstar,sigAstar
	//else
		//return sigU,sigA
	u = gsl_rng_uniform(self->rand);
	if(self->prob<u){
		gsl_vector_memcpy(self->sigU,self->sigUstar);
		gsl_vector_memcpy(self->sigA,self->sigAstar);
	}

	return(0);
GMERRH("NGPmcmc_drawSigs",1);
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
#undef DDOT_VVOUTER

/*************************************************
* Class method for doing the operation once. It  *
* basically functions like what a function looks *
* like for a normal person with weird notation.  *
* If the operation is going to be performed more *
* than once do not use this as you will use more *
* time for memory allocation 					 *
*************************************************/

int NGPmcmc_NGPmcmc(gsl_matrix *y,gsl_vector *tobs,int Niter,gsl_vector *sigU,
					gsl_vector *sigA,double sigEps,double sigMu,double sigAlph,
					double a,double b,marray3d *th,gsl_matrix *sig){
	NGPmcmc *self=NULL;
	if((self = NGPmcmc_New())==NULL)	GMERR(-1);
	if(NGPmcmc_construct(self,y,tobs,Niter,sigU,sigA,sigEps,sigMu,sigAlph,a,b,th,sig))			GMERR(-11);
	if(NGPmcmc_operations(self))		GMERR(-21);
	if(NGPmcmc_free(self))				GMERR(-31);
	return(0);
GMERRH("NGPmcmc_NGPmcmc",1);
}

/**********************************************
* Any parameter modification methods go here. *
* Make sure when these are called that the	  *
* outputs are zeroed out					  *
**********************************************/

int NGPmcmc_updatePrior(NGPmcmc *self,char parameter,gsl_vector *param){
	if((self==NULL)||(param==NULL))		GMERR(-1);
	if(_NGPmcmc_zeroOutputs(self))		GMERR(-11);
	switch(parameter){
		case 'a':
			self->a = gsl_vector_get(param,0);
			printf("Reset a=%f\n",gsl_vector_get(param,0));
			break;
		case 'b':
			self->b = gsl_vector_get(param,0);
			printf("Reset b=%f\n",gsl_vector_get(param,0));
			break;
		case 'm':
			self->sigMu = gsl_vector_get(param,0);
			printf("Reset sigMu=%f\n",gsl_vector_get(param,0));
			break;
		case 'l':
			self->sigAlph = gsl_vector_get(param,0);
			printf("Reset sigAlph=%f\n",gsl_vector_get(param,0));
			break;
		default:
			GMERR(-21);

	}
	return(0);
GMERRH("NGPmcmc_updatePrior",1);
}


/************************************************
* This section writes the outputs of the method *
* to a textfile	using NGPmcmc_writeOutputs		* 
************************************************/

int NGPmcmc_writeTheta(NGPmcmc *self,char *name){
	FILE *fp;
	if(self==NULL)					GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)					GMERR(-3);
	if(marray3d_write(fp,self->th))	GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPmcmc_writeTheta",1);
}

int NGPmcmc_writeSig(NGPmcmc *self,char *name){
	FILE *fp;
	if(self==NULL)					GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)					GMERR(-3);
	if(gsl_matrix_fprintf(fp,self->sig,"%0.6f"))	GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPmcmc_writeSig",1);
}

int NGPmcmc_writeOutputs(NGPmcmc *self,char *baseName){
	//Here we define all the names for each output
	char strTheta[16] = "_Theta.txt";
	char outTheta[1024];
	char strSig[16] = "_Sig.txt";
	char outSig[1024];
	sprintf(outTheta,"%s%s",baseName,strTheta);
	sprintf(outSig,"%s%s",baseName,strSig);
	
	printf("%s\n",outTheta);
	printf("%s\n",outSig);
	//Here we actually write out all the data
	if(self==NULL)							GMERR(-1);
	if(NGPmcmc_writeTheta(self,outTheta))	GMERR(-11);
	if(NGPmcmc_writeSig(self,outSig))		GMERR(-21);
	return(0);
GMERRH("NGPmcmc_writeOutputs",1);
}


#undef NGPmcmc_MIN
