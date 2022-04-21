#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include "gmlib.h"
#include "matrix_utils.h"
#include "mvnorm.h"
#include "mdarray.h"
#include "KalmanFilter2.h"
#include "fastStateSmoother2.h"
#include "SSsimulate2.h"
#include "NGPtools2.h"
#include "NGPsample2.h"
#include "NGPlaf2.h"

//Function declarations
static int _NGPlaf2_checkInputs(gsl_vector *tobs,gsl_matrix *y,int NK,int NL,int Niter,marray3d *Ksi_out,marray3d *Psi_out,marray3d *yhat_out,marray3d *mu_out,marray4d *Sigma_out,gsl_vector *sigPrior,gsl_vector *epsPrior,gsl_vector *ksiPrior,gsl_vector *APrior,gsl_vector *psiPrior,gsl_vector *BPrior,gsl_vector *aaPrior);
int NGPlaf2_sample_Ksi(NGPlaf2 *self);
int NGPlaf2_sample_sigKsi(NGPlaf2 *self);
int NGPlaf2_sample_sigA(NGPlaf2 *self);
int NGPlaf2_sample_Sigma0(NGPlaf2 *self);
int NGPlaf2_sample_Psi(NGPlaf2 *self);
int NGPlaf2_sample_sigPsi(NGPlaf2 *self);
int NGPlaf2_sample_sigB(NGPlaf2 *self);
int NGPlaf2_sample_eta(NGPlaf2 *self);
int NGPlaf2_sample_Theta(NGPlaf2 *self);
int NGPlaf2_sample_Phi(NGPlaf2 *self);
int NGPlaf2_sample_Tau(NGPlaf2 *self);

NGPlaf2 * NGPlaf2_New(){
	NGPlaf2 *self;
	self = (NGPlaf2 *)GM_Malloc(sizeof(NGPlaf2));
	if(self==NULL)								GMERR(-1);
	if(NGPlaf2_init(self))						GMERR(-21);
	return(self);
GMERRH("NGPlaf2_New",NULL);
}

#define NNN(A) self->A=NULL

int NGPlaf2_init(NGPlaf2 *self){
	if(self==NULL)								GMERR(-1);
	memset(self,'\0',sizeof(NGPlaf2));

	self->Np=0;
	self->Nt=0;
	self->NK=0;
	self->NL=0;
	self->Nst=0; self->Nm=0;

	//prior parameters
	self->sigMu = 0.0;
	self->sigAlph = 0.0;
	self->aEps = 0.0;
	self->bEps = 0.0;
	self->aKsi = 0.0;
	self->bKsi = 0.0;
	self->aA = 0.0;
	self->bA = 0.0;
	self->aPsi = 0.0;
	self->bPsi = 0.0;
	self->aB = 0.0;
	self->bB = 0.0;
	self->a1 = 0.0;
	self->a2 = 0.0;

	self->Niter = 0;
	self->burnin = 0;
	self->Nthin = 0;

	/*********
	* Inputs *
	*********/
	NNN(tobs);
	NNN(y);
	NNN(delta);


	/**********
	* Outputs *
	**********/

	NNN(Ksi_out);

	/*****************
	* Main variables *
	*****************/
	NNN(eta);
	NNN(Theta);
	NNN(Phi);
	NNN(theta);
	NNN(tau);
	NNN(Sigma0);
	NNN(sigKsi);
	NNN(sigA);
	NNN(sigPsi);
	NNN(sigB);
	NNN(DKsi);
	NNN(DPsi);
	NNN(B);
	NNN(A);

	/*************
	* sample_Ksi *
	*************/
	//sampler
	NNN(ZDummy_Ksi);
	NNN(HDummy_Ksi);
	NNN(etat_Ksi);
	NNN(ketatTheta_Ksi);
	NNN(Z_Ksi);
	NNN(H_Ksi);
	NNN(T_Ksi);
	NNN(R_Ksi);
	NNN(Q_Ksi);
	NNN(P0_Ksi);
	NNN(a0_Ksi);
	NNN(KsiSampler);

	/****************
	* sample_sigKsi *
	****************/
	NNN(RKsi_sigKsi);

	/**************
	* sample_sigA *
	**************/
	NNN(RA_sigA);

	/****************
	* sample_Sigma0 *
	****************/
	//vectors
	NNN(etat_Sigma0);
	NNN(Ksit_Sigma0);
	NNN(yhatt_Sigma0);
	NNN(SS_Sigma0);

	//matrices
	NNN(kron_Sigma0);
	NNN(yhat_Sigma0);

	/*************
	* sample_Psi *
	*************/
	NNN(ZDummy_Psi);
	NNN(HDummy_Psi);
	NNN(Ksit_Psi);
	NNN(KsitM_Psi);
	NNN(Z_t_Psi);
	NNN(Z_tT_Psi);
	NNN(ZZT_Psi);

	//used for the sampler
	NNN(Z_Psi);
	NNN(H_Psi);
	NNN(T_Psi);
	NNN(R_Psi);
	NNN(Q_Psi);
	NNN(alphaSamples_Psi);
	NNN(P0_Psi);
	NNN(a0_Psi);
	NNN(PsiSampler);


	/****************
	* sample_sigPsi *
	****************/
	NNN(RPsi_sigPsi);

	/**************
	* sample_sigB *
	**************/
	NNN(RB_sigB);

	/*************
	* sample_eta *
	*************/
	//Vectors
	NNN(Ksit_Eta);
	NNN(Psit_Eta);
	NNN(yhat_Eta);
	NNN(yt_Eta);
	NNN(mu_Eta);
	NNN(etat_Eta);

	//Matrices
	NNN(KsitM_Eta);
	NNN(Z_t_Eta);
	NNN(Z_tT_Eta);
	NNN(Lambda_Eta);
	NNN(Sigma_Eta);
	NNN(SigZ_Eta);
	NNN(eye_Eta);

	/***************
 	* sample_Theta *
	***************/
	//vectors
	NNN(Ksit_Theta);
	NNN(etat_Theta);
	NNN(eta_tildet_Theta);
	NNN(Phit_Theta);
	NNN(yp_Theta);
	NNN(etay_Theta);
	NNN(mu_Theta);
	NNN(Thetat_Theta);

	//matrices
	NNN(KsitM_Theta);
	NNN(eta_tilde_Theta);
	NNN(eta_tildeT_Theta);
	NNN(etaeta_Theta);
	NNN(Siginv_Theta);
	NNN(SigTheta_Theta);
	NNN(tol_Theta);
	NNN(Phitau_Theta);

	/*************
	* sample_Phi *
	*************/
	//None

	/*************
	* sample_tau *
	*************/

	NNN(Theta2_Tau);
	NNN(tau_minus_Tau);
	NNN(sum_Tau);
	
	/**************************
	* random number generator *
	**************************/
	NNN(rand);
	
	return(0);
GMERRH("NGPlaf2_init",1);
}

#undef NNN




static int _NGPlaf2_checkInputs(gsl_vector *tobs,gsl_matrix *y,int NK,int NL,int Niter,marray3d *Ksi_out,
								marray3d *Psi_out,marray3d *yhat_out,marray3d *mu_out,
								marray4d *Sigma_out,gsl_vector *sigPrior,gsl_vector *epsPrior,
								gsl_vector *ksiPrior,gsl_vector *APrior,gsl_vector *psiPrior,
								gsl_vector *BPrior,gsl_vector *aaPrior){
	int Nt,Np,Nsamps;

	if((tobs==NULL)||(y==NULL)||(Ksi_out==NULL)||(Psi_out==NULL)||
		(yhat_out==NULL)||(mu_out==NULL)||(Sigma_out==NULL)||(sigPrior==NULL)||
		(epsPrior==NULL)||(ksiPrior==NULL)||(APrior==NULL)||
		(psiPrior==NULL)||(BPrior==NULL)||(aaPrior==NULL))	GMERR(-1);
	Nsamps = Niter/5;

	//All these vectors have the prior parameters for 
	//the gamma distributions with 2 elements each
	if(sigPrior->size!=2)					GMERR(-11);
	if(epsPrior->size!=2)					GMERR(-21);
	if(ksiPrior->size!=2)					GMERR(-31);
	if(APrior->size!=2)						GMERR(-41);
	if(psiPrior->size!=2)					GMERR(-41);
	if(BPrior->size!=2)						GMERR(-51);
	if(aaPrior->size!=2)					GMERR(-52);
	
	//All the priors should be positive
	if(positiveVector(sigPrior))			GMERR(-61);
	if(positiveVector(epsPrior))			GMERR(-71);
	if(positiveVector(ksiPrior))			GMERR(-81);
	if(positiveVector(psiPrior))			GMERR(-91);
	if(positiveVector(APrior))				GMERR(-101);
	if(positiveVector(BPrior))				GMERR(-111);
	if(positiveVector(aaPrior))				GMERR(-121);

	//Check sizes
	Nt = y->size1;
	Np = y->size2;

	if(tobs->size!=Nt)						GMERR(-131);

	//Make sure that the output dimensions are fine
	if(marray3d_checkSizes(Ksi_out,Nsamps,Nt,NK*NL))		GMERR(-141);
	if(marray3d_checkSizes(Psi_out,Nsamps,Nt,NK))			GMERR(-151);
	if(marray3d_checkSizes(mu_out,Nsamps,Nt,Np))			GMERR(-161);
	if(marray3d_checkSizes(yhat_out,Nsamps,Nt,Np))			GMERR(-171);
	if(marray4d_checkSizes(Sigma_out,Nsamps,Nt,Np,Np))		GMERR(-181);

	return(0);
GMERRH("_NGPlaf2_checkInputs",1);
}

#define M_ALLOC(M,SIZE1,SIZE2) self->M = gsl_matrix_alloc(SIZE1,SIZE2)
#define V_ALLOC(V,SIZE) self->V = gsl_vector_alloc(SIZE)
#define A_ALLOC3(A,SIZE1,SIZE2,SIZE3) self->A = marray3d_alloc(SIZE1,SIZE2,SIZE3)
#define A_ALLOC4(A,SIZE1,SIZE2,SIZE3,SIZE4) self->A = marray4d_alloc(SIZE1,SIZE2,SIZE3)

int NGPlaf2_construct(NGPlaf2 *self,gsl_vector *tobs,gsl_matrix *y,int NK,int NL,int Niter,marray3d *Ksi_out,marray3d *Psi_out,marray3d *yhat_out,marray3d *mu_out,marray4d *Sigma_out,gsl_vector *sigPrior,gsl_vector *epsPrior,gsl_vector *ksiPrior,gsl_vector *APrior,gsl_vector *psiPrior,gsl_vector *BPrior,gsl_vector *aaPrior){
	int Np,Nt,Nst,Nm;
	int i;
	double temp;
	if(self==NULL)								GMERR(-1);
	if(_NGPlaf2_checkInputs(tobs,y,NK,NL,Niter,Ksi_out,Psi_out,yhat_out,
						mu_out,Sigma_out,sigPrior,epsPrior,
						ksiPrior,APrior,psiPrior,BPrior,aaPrior))			GMERR(-11);

	self->NK = NK;
	self->NL = NL;
	self->Nt = y->size1;
	self->Np = y->size2;
	self->Nst = NK*NL;
	Nst = self->Nst;
	self->Nm = 3*Nst;
	Nm = self->Nm;
	self->Niter=Niter;
	self->burnin = 5000;
	self->Nthin = 5;

	Np = self->Np;
	Nt = self->Nt;

	/*********
	* Inputs *
	*********/
	self->y = y;
	self->tobs = tobs;
	V_ALLOC(delta,Nt);

	//Taking the differences. Not sure what to do with 
	//the last one but I'm not sure if the algorithm needs it
	for(i=0;i<Nt-1;i++){
		temp = gsl_vector_get(tobs,i+1)-gsl_vector_get(tobs,i);
		gsl_vector_set(self->delta,i,temp);
	}
	gsl_vector_set(self->delta,Nt-1,temp);



	/**************************
	* random number generator *
	**************************/
	self->rand = gsl_rng_alloc(gsl_rng_ranlxs0);

	//prior parameters
	self->sigMu = gsl_vector_get(sigPrior,0);
	self->sigAlph = gsl_vector_get(sigPrior,1);
	self->aEps = gsl_vector_get(epsPrior,0);
	self->bEps = gsl_vector_get(epsPrior,1);
	self->aKsi = gsl_vector_get(ksiPrior,0);
	self->bKsi = gsl_vector_get(ksiPrior,1);
	self->aA = gsl_vector_get(APrior,0);
	self->bA = gsl_vector_get(APrior,1);
	self->aPsi = gsl_vector_get(psiPrior,0);
	self->bPsi = gsl_vector_get(psiPrior,1);
	self->aB = gsl_vector_get(BPrior,0);
	self->bB = gsl_vector_get(BPrior,1);
	self->a1 = gsl_vector_get(aaPrior,0);
	self->a2 = gsl_vector_get(aaPrior,1);

	/********************
	* overall variables *
	********************/
	//eta = 0.1 * randn(NK,1,Nt)
	//Theta = 0.1 * randn(Np,NL)
	//Phi = rand(Gamma(3/2,2/3),Np,NL)
	//theta = rand(Gamma(10,1),NL)
	//tau = cumprod(theta)
	//Sigma0 = diagm(rand(InverseGamma(1,0.1),Np))
	//sigKsi = 50. * ones(NL *NK)
	//sigA  = 50. * ones(NL * NK)
	//sigPsi = 50. * ones(NK)
	//sigB = 50. * ones(NK)
	self->sigKsi = gsl_vector_set_ones(NL*NK);
	self->sigA = gsl_vector_set_ones(NL*NK);
	self->sigPsi = gsl_vector_set_ones(NK);
	self->sigB = gsl_vector_set_ones(NK);
	gsl_vector_scale(self->sigA,50.0);
	gsl_vector_scale(self->sigPsi,50.0);
	gsl_vector_scale(self->sigB,50.0);
	gsl_vector_scale(self->sigKsi,50.0);
	M_ALLOC(eta,Nt,NK);
	M_ALLOC(Theta,Np,NL);
	M_ALLOC(Phi,Np,NL);
	V_ALLOC(theta,NL);
	V_ALLOC(tau,NL);
	M_ALLOC(Sigma0,Np,Np);
	M_ALLOC(DKsi,Nt,Nst);
	M_ALLOC(A,Nt,Nst);
	M_ALLOC(DPsi,Nt,NK);
	M_ALLOC(B,Nt,NK);
	M_ALLOC(Psi,Nt,NK);
	M_ALLOC(Ksi,Nt,NK*NL);
	if(randnM(self->eta,self->rand,0,1))				GMERR(-21);
	if(randnM(self->Theta,self->rand,0,1))				GMERR(-31);
	if(randgM(self->Phi,self->rand,1.5,.667))			GMERR(-41);
	if(randgV(self->theta,self->rand,10.0,1.0))			GMERR(-51);
	if(cumProd(self->theta,self->tau))					GMERR(-61);
	if(randDiagIG(self->Sigma0,self->rand,1.0,0.1))		GMERR(-71);
	if(gsl_matrix_scale(self->eta,0.1))					GMERR(-81);
	if(gsl_matrix_scale(self->Theta,0.1))				GMERR(-91);

	/*************
	* sample_Ksi *
	*************/
	//vectors
	V_ALLOC(etat_Ksi,NK);
	V_ALLOC(a0_Ksi,Nm);
	V_ALLOC(work_Ksi,Nt);

	//matrices
	M_ALLOC(P0_Ksi,Nm,Nm);
	M_ALLOC(ketatTheta_Ksi,Np,NK*NL);

	//arrays
	A_ALLOC3(ZDummy_Ksi,Nt,Np,Nm);
	A_ALLOC3(HDummy_Ksi,Nt,Np,Np);
	A_ALLOC3(Z_Ksi,Nt,Np,Nm);
	A_ALLOC3(H_Ksi,1,Np,Np);
	A_ALLOC3(T_Ksi,Nt,Nm,Nm);
	A_ALLOC3(R_Ksi,Nt,Nm,2*Nst);
	A_ALLOC3(Q_Ksi,Nt,2*Nst,2*Nst);
	A_ALLOC3(alphaSamples_Ksi,1,Nt,Nm);

	//sampler
	self->KsiSampler = _NGPsample_New();
	if(_NGPsample_construct(self->KsiSampler,self->y,1,self->a0_Ksi,
							self->P0_Ksi,self->Z_Ksi,self->H_Ksi,
							self->T_Ksi,self->R_Ksi,self->Q_Ksi,
							self->alphaSamples_Ksi))	GMERR(-101);

	/****************
	* sample_sigKsi *
	****************/
	V_ALLOC(RKsi_sigKsi,NK*NL);

	/**************
	* sample_sigA *
	**************/
	V_ALLOC(RA_sigA,NK*NL);

	/****************
	* sample_Sigma0 *
	****************/
	//vectors
	V_ALLOC(etat_Sigma0,NK);
	V_ALLOC(Ksit_Sigma0,NL*NK);
	V_ALLOC(yhatt_Sigma0,Np);
	V_ALLOC(SS_Sigma0,Np);

	//matrices
	M_ALLOC(kron_Sigma0,Np,NL*NK);
	M_ALLOC(yhat_Sigma0,Nt,Np);

	/*************
	* sample_Psi *
	*************/
	//vectors
	V_ALLOC(Ksit_Psi,NL*NK);
	V_ALLOC(a0_Psi,3*NK);
	V_ALLOC(work_Psi,Nt);

	//matrices
	M_ALLOC(Z_t_Psi,Np,NK);
	M_ALLOC(Z_tT_Psi,NK,Np);
	M_ALLOC(ZZT_Psi,Np,Np);
	M_ALLOC(P0_Psi,3*NK,3*NK);
	M_ALLOC(KsitM_Psi,NL,NK);

	//arrays
	A_ALLOC3(ZDummy_Psi,Nt,Np,3*NK);
	A_ALLOC3(HDummy_Psi,Nt,Np,Np);
	A_ALLOC3(Z_Psi,Nt,Np,3*NK);
	A_ALLOC3(H_Psi,Nt,Np,Np);
	A_ALLOC3(T_Psi,Nt,3*NK,3*NK);
	A_ALLOC3(R_Psi,Nt,3*NK,2*NK);
	A_ALLOC3(Q_Psi,Nt,2*NK,2*NK);
	A_ALLOC3(alphaSamples_Psi,1,Nt,3*NK);

	//sampler
	self->PsiSampler = _NGPsample_New();
	if(_NGPsample_construct(self->PsiSampler,self->y,1,self->a0_Psi,
							self->P0_Psi,self->Z_Psi,self->H_Psi,
							self->T_Psi,self->R_Psi,self->Q_Psi,
							self->alphaSamples_Psi))	GMERR(-111);

	/****************
	* sample_sigPsi *
	****************/
	V_ALLOC(RPsi_sigPsi,NK);

	/**************
	* sample_sigB *
	**************/
	V_ALLOC(RB_sigB,NK);

	/*************
	* sample_eta *
	*************/
	//vectors
	V_ALLOC(etat_Eta,NK);
	V_ALLOC(Ksit_Eta,NL*NK);
	V_ALLOC(mu_Eta,NK);
	V_ALLOC(Psit_Eta,NK);
	V_ALLOC(yhat_Eta,Np);
	V_ALLOC(yt_Eta,Np);

	//matrices
	M_ALLOC(eye_Eta,NK,NK);
	gsl_matrix_set_identity(self->eye_Eta);
	M_ALLOC(KsitM_Eta,NL,NK);
	M_ALLOC(Lambda_Eta,NK,NK);
	M_ALLOC(Sigma_Eta,NK,NK);
	M_ALLOC(SigZ_Eta,NK,Np);
	M_ALLOC(Z_t_Eta,Np,NK);
	M_ALLOC(Z_tT_Eta,NK,Np);

	/***************
 	* sample_Theta *
	***************/
	//vectors
	V_ALLOC(eta_tildet_Theta,NK);
	V_ALLOC(etat_Theta,NK);
	V_ALLOC(etay_Theta,NL);
	V_ALLOC(Ksit_Theta,NL*NK);
	V_ALLOC(mu_Theta,NL);
	V_ALLOC(Phit_Theta,NL);
	V_ALLOC(Thetat_Theta,NL);
	V_ALLOC(yp_Theta,Nt);

	//matrices
	M_ALLOC(eta_tilde_Theta,Nt,NL);
	M_ALLOC(eta_tildeT_Theta,NL,Nt);
	M_ALLOC(etaeta_Theta,NL,NL);
	M_ALLOC(KsitM_Theta,NL,NK);
	M_ALLOC(Siginv_Theta,NL,NL);
	M_ALLOC(SigTheta_Theta,NL,NL);
	M_ALLOC(Phitau_Theta,NL,NL);
	gsl_matrix_set_zero(self->Phitau_Theta);

	//tol_Theta not updated
	M_ALLOC(tol_Theta,NL,NL);
	gsl_matrix_set_identity(self->tol_Theta);
	gsl_matrix_scale(self->tol_Theta,.00000001);

	/*************
	* sample_Phi *
	*************/
	//NONE

	/*************
	* sample_tau *
	*************/
	//vectors
	V_ALLOC(sum_Tau,NL);
	V_ALLOC(tau_minus_Tau,NL);

	//matrices
	M_ALLOC(Theta2_Tau,Np,NL);

	/*************
	* Operations *
	*************/
	V_ALLOC(eta_t_O,NK);
	V_ALLOC(Ksi_t_O,NL*NK);
	V_ALLOC(mu_t_O,Np);
	V_ALLOC(yhat_t_O,Np);
	V_ALLOC(Psi_t_O,NK);

	M_ALLOC(Ksi_this_O,NL,NK);
	M_ALLOC(ThetaKsi_O,Np,NK);
	M_ALLOC(ThetaKsiT_O,NK,Np);
	M_ALLOC(TKKT_O,Np,Np);

	/**********
	* outputs *
	**********/
	self->Ksi_out = Ksi_out;
	self->Psi_out = Psi_out;
	self->yhat_out = yhat_out;
	self->mu_out = mu_out;
	self->Sigma_out = Sigma_out;

	return(0);
GMERRH("NGPlaf2_construct",1);
}

int NGPlaf2_reset(NGPlaf2 *self){
	int NL,NK;
	if(self==NULL)								GMERR(-1);
	NL = self->NL;NK = self->NK;
	self->sigKsi = gsl_vector_set_ones(NL*NK);
	self->sigA = gsl_vector_set_ones(NL*NK);
	self->sigPsi = gsl_vector_set_ones(NK);
	self->sigB = gsl_vector_set_ones(NK);
	gsl_vector_scale(self->sigA,50.0);
	gsl_vector_scale(self->sigPsi,50.0);
	gsl_vector_scale(self->sigB,50.0);
	gsl_vector_scale(self->sigKsi,50.0);
	if(randnM(self->eta,self->rand,0,1))				GMERR(-21);
	if(randnM(self->Theta,self->rand,0,1))				GMERR(-31);
	if(randgM(self->Phi,self->rand,1.5,.667))			GMERR(-41);
	if(randgV(self->theta,self->rand,10.0,1.0))			GMERR(-51);
	if(cumProd(self->theta,self->tau))					GMERR(-61);
	if(randDiagIG(self->Sigma0,self->rand,1.0,0.1))		GMERR(-71);
	if(gsl_matrix_scale(self->eta,0.1))					GMERR(-81);
	if(gsl_matrix_scale(self->Theta,0.1))				GMERR(-91);
	marray3d_set_zero(self->Ksi_out);
	marray3d_set_zero(self->Psi_out);
	marray3d_set_zero(self->yhat_out);
	marray3d_set_zero(self->mu_out);
	marray4d_set_zero(self->Sigma_out);
	return(0);
GMERRH("NGPlaf2_resetSigma",1);
}

int NGPlaf2_modifyPrior(NGPlaf2 *self,char parameter,gsl_vector *param){
	if(self==NULL)								GMERR(-1);
	switch(parameter){
		case 'a':
			self->a1 = gsl_vector_get(param,0);
			self->a2 = gsl_vector_get(param,1);
			printf("Reset a1=%f,a2=%f\n",gsl_vector_get(param,0),gsl_vector_get(param,1));
			break;
		case 'A':
			self->aA = gsl_vector_get(param,0);
			self->bA = gsl_vector_get(param,1);
			printf("Reset aA=%f,bA=%f\n",gsl_vector_get(param,0),gsl_vector_get(param,1));
			break;
		case 'B':
			self->aB = gsl_vector_get(param,0);
			self->bB = gsl_vector_get(param,1);
			printf("Reset aB=%f,aB=%f\n",gsl_vector_get(param,0),gsl_vector_get(param,1));
			break;
		case 'S':
			self->sigMu = gsl_vector_get(param,0);
			self->sigAlph = gsl_vector_get(param,1);
			printf("Reset sigMu=%f,sigAlph=%f\n",gsl_vector_get(param,0),gsl_vector_get(param,1));
			break;
		case 'K':
			self->aKsi = gsl_vector_get(param,0);
			self->bKsi = gsl_vector_get(param,1);
			printf("Reset aKsi=%f,bKsi=%f\n",gsl_vector_get(param,0),gsl_vector_get(param,1));
			break;
		case 'P':
			self->aPsi = gsl_vector_get(param,0);
			self->bPsi = gsl_vector_get(param,1);
			printf("Reset aPsi=%f,bPsi=%f\n",gsl_vector_get(param,0),gsl_vector_get(param,1));
			break;
		case 'E':
			self->aEps = gsl_vector_get(param,0);
			self->bEps = gsl_vector_get(param,1);
			printf("Reset aEps=%f,bEps=%f\n",gsl_vector_get(param,0),gsl_vector_get(param,1));
			break;
		default:
			GMERR(-11);
	}
	return(0);
GMERRH("NGPlaf2_modifyPrior",1);
}

#undef V_ALLOC
#undef M_ALLOC
#undef A_ALLOC3
#undef A_ALLOC4

int NGPlaf2_free(NGPlaf2 *s){
	if(s==NULL)									GMERR(-1);
	if(&GM_FreeGSLVector==NULL)					GMERR(-2);
	if(&GM_FreeGSLMatrix==NULL)					GMERR(-3);

	/**********
	* overall *
	**********/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->sigKsi,&s->sigA,&s->sigB,
				&s->sigPsi,&s->tau,&s->theta,NULL))		GMERR(-1011);


	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
				&s->Theta,&s->Phi,&s->eta,
				&s->A,&s->B,&s->DKsi,&s->DPsi,
				&s->Sigma0,&s->Ksi,
				&s->Psi,NULL))							GMERR(-1021);

	/***************
	* sample_Theta *
	***************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->Ksit_Theta,&s->etat_Theta,
				&s->eta_tildet_Theta,&s->etay_Theta,
				&s->Phit_Theta,&s->yp_Theta,
				&s->mu_Theta,&s->Thetat_Theta,NULL))	GMERR(-311);
	
	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
				&s->KsitM_Theta,&s->eta_tilde_Theta,
				&s->eta_tildeT_Theta,&s->etaeta_Theta,
				&s->Siginv_Theta,&s->tol_Theta,
				&s->Phitau_Theta,
				&s->SigTheta_Theta,NULL))				GMERR(-321);

	/*************
	* sample_Eta *
	*************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->Ksit_Eta,&s->Psit_Eta,
				&s->yhat_Eta,&s->yt_Eta,
				&s->mu_Eta,&s->etat_Eta,NULL))			GMERR(-511);	
	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
				&s->KsitM_Eta,&s->Z_t_Eta,
				&s->Z_tT_Eta,&s->Lambda_Eta,
				&s->Sigma_Eta,&s->SigZ_Eta,
				&s->eye_Eta,NULL))						GMERR(-521);

	/*************
	* sample_Phi *
	*************/
	//None

	/*************
	* sample_Ksi *
	*************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->work_Ksi,&s->a0_Ksi,&s->etat_Ksi,NULL))	GMERR(-1311);

	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
				&s->P0_Ksi,&s->ketatTheta_Ksi,NULL))	GMERR(-1321);

	marray3d_free(s->ZDummy_Ksi);
	marray3d_free(s->HDummy_Ksi);
	marray3d_free(s->Z_Ksi);
	marray3d_free(s->H_Ksi);
	marray3d_free(s->Q_Ksi);
	marray3d_free(s->T_Ksi);
	marray3d_free(s->R_Ksi);
	marray3d_free(s->alphaSamples_Ksi);
	
	if(_NGPsample_free(s->KsiSampler))					GMERR(-1341);

	/****************
	* sample_sigKsi *
	****************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->RKsi_sigKsi,NULL))					GMERR(-811);

	/**************
	* sample_sigA *
	**************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->RA_sigA,NULL))						GMERR(-711);
	
	/*************
	* sample_Psi *
	*************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",&s->work_Psi,
				&s->Ksit_Psi,&s->a0_Psi,NULL))			GMERR(-1211);

	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
				&s->KsitM_Psi,&s->Z_t_Psi,
				&s->Z_tT_Psi,&s->ZZT_Psi,
				&s->P0_Psi,NULL))						GMERR(-1221);

	marray3d_free(s->ZDummy_Psi);
	marray3d_free(s->HDummy_Psi);
	marray3d_free(s->Z_Psi);
	marray3d_free(s->H_Psi);
	marray3d_free(s->T_Psi);
	marray3d_free(s->R_Psi);
	marray3d_free(s->Q_Psi);
	marray3d_free(s->alphaSamples_Psi);

	
	if(_NGPsample_free(s->PsiSampler))					GMERR(-1241);

	/****************
	* sample_sigPsi *
	****************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->RPsi_sigPsi,NULL))					GMERR(-911);
	
	/**************
	* sample_sigB *
	**************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->RB_sigB,NULL))						GMERR(-611);

	/****************
	* sample_Sigma0 *
	****************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->etat_Sigma0,&s->Ksit_Sigma0,
				&s->yhatt_Sigma0,
				&s->SS_Sigma0,NULL))					GMERR(-211);
	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
				&s->kron_Sigma0,&s->yhat_Sigma0,NULL))	GMERR(-221);
	
	/*************
	* sample_Tau *
	*************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->sum_Tau,&s->tau_minus_Tau,NULL))	GMERR(-411);
	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
				&s->Theta2_Tau,NULL))					GMERR(-421);

	/*************
	* Operations *
	*************/
	if(GM_FreezeListWithMethod(&GM_FreeGSLVector,"x",
				&s->eta_t_O,&s->Ksi_t_O,
				&s->mu_t_O,&s->Psi_t_O,
				&s->yhat_t_O,NULL))						GMERR(-1411);
	if(GM_FreezeListWithMethod(&GM_FreeGSLMatrix,"x",
				&s->Ksi_this_O,
				&s->ThetaKsi_O,
				&s->ThetaKsiT_O,
				&s->TKKT_O,NULL))						GMERR(-1421);

	return(0);
GMERRH("NGPlaf2_free",1);
}

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

int NGPlaf2_operations(NGPlaf2 *self){
	int t,it,burn,iter,Nt,NL,NK,thin;
	int count = 0;
	if(self==NULL)								GMERR(-1);
	burn = self->burnin;
	iter = self->Niter;
	Nt = self->Nt;
	NL = self->NL;
	NK = self->NK;
	thin = self->Nthin;
	
	for(it=0;it<burn;it++){
		if((it%100)==0)	printf("Burn %d\n",it);
		if(NGPlaf2_sample_Ksi(self))					GMERR(-11);
		if(NGPlaf2_sample_sigKsi(self))					GMERR(-21);
		if(NGPlaf2_sample_sigA(self))					GMERR(-31);
		if(NGPlaf2_sample_Psi(self))					GMERR(-41);
		if(NGPlaf2_sample_sigPsi(self))					GMERR(-51);
		if(NGPlaf2_sample_sigB(self))					GMERR(-61);
		if(NGPlaf2_sample_eta(self))					GMERR(-71);
		if(NGPlaf2_sample_Sigma0(self))					GMERR(-81);
		if(NGPlaf2_sample_Theta(self))					GMERR(-91);
		if(NGPlaf2_sample_Phi(self))					GMERR(-101);
		if(NGPlaf2_sample_Tau(self))					GMERR(-111);
	}
	for(it=0;it<iter;it++){
		if((it%100)==0)	printf("Iteration %d\n",it);
		if(NGPlaf2_sample_Ksi(self))					GMERR(-211);
		if(NGPlaf2_sample_sigKsi(self))					GMERR(-221);
		if(NGPlaf2_sample_sigA(self))					GMERR(-231);
		if(NGPlaf2_sample_Psi(self))					GMERR(-241);
		if(NGPlaf2_sample_sigPsi(self))					GMERR(-251);
		if(NGPlaf2_sample_sigB(self))					GMERR(-261);
		if(NGPlaf2_sample_eta(self))					GMERR(-271);
		if(NGPlaf2_sample_Sigma0(self))					GMERR(-281);
		if(NGPlaf2_sample_Theta(self))					GMERR(-291);
		if(NGPlaf2_sample_Phi(self))					GMERR(-301);
		if(NGPlaf2_sample_Tau(self))					GMERR(-311);
		
		//Here is where we log our results
		if((it%thin)==0){
			//Ksi_mat = Ksi
			if(marray3d_set_X(self->Ksi_out,self->Ksi,count))	GMERR(-312);

			//Psi_mat = Psi
			if(marray3d_set_X(self->Psi_out,self->Psi,count))	GMERR(-313);

			for(t=0;t<Nt;t++){
				//Ksi_this = reshape(Ksi[:,t],NL,NK)
				if(M_GET_ROW(Ksi,Ksi_t_O,t))							GMERR(-321);
				if(reshape(self->Ksi_t_O,NL,NK,self->Ksi_this_O))		GMERR(-331);

				//yhat[:,t,it] = Theta * Ksi_this * slice(eta,:,1,t)
				if(DDOT_MM(Theta,Ksi_this_O,ThetaKsi_O))				GMERR(-341);
				if(M_GET_ROW(eta,eta_t_O,t))							GMERR(-351);
				if(DDOT_MV(ThetaKsi_O,eta_t_O,yhat_t_O))				GMERR(-361);
				if(marray3d_set_pencilZ(self->yhat_out,self->yhat_t_O,count,t))GMERR(-371);

				//mu[:,t,it] = Theta * Ksi_this * slice(Psi,:,t)
				if(M_GET_ROW(Psi,Psi_t_O,t))							GMERR(-381);
				if(DDOT_MV(ThetaKsi_O,Psi_t_O,mu_t_O))					GMERR(-391);
				if(marray3d_set_pencilZ(self->mu_out,self->mu_t_O,count,t))GMERR(-401);

				//Sigma[:,:,t,it] = Theta * Ksi_this * Ksi_this' * Theta' + Sigma0
				if(gTranspose(self->ThetaKsi_O,self->ThetaKsiT_O))		GMERR(-411);
				if(DDOT_MM(ThetaKsi_O,ThetaKsiT_O,TKKT_O))				GMERR(-421);
				if(M_ADD(TKKT_O,Sigma0))								GMERR(-431);
				if(marray4d_set_ZW(self->Sigma_out,self->TKKT_O,count,t))	GMERR(-441);
			}
			count = count + 1;
		}
	}
	
	return(0);
GMERRH("NGPlaf2_operations",1);
}

int NGPlaf2_sample_Ksi(NGPlaf2 *self){
	int t,Nt,j,Nst,k,Np;
	Nt = self->Nt;
	double temp;
	Nst = self->Nst;
	Np = self->Np;
	double zeroEps = 0;
	
	//Z = zeros(Float64,Np,Nm,Nt)
	marray3d_set_zero(self->Z_Ksi);
	for(t=0;t<Nt;t++){
		//Z[:,1:Nst,t] = kron(slice(eta,:,:,t)',Theta) 
		if(M_GET_ROW(eta,etat_Ksi,t))									GMERR(-11);
		if(kron_vm(self->etat_Ksi,self->Theta,self->ketatTheta_Ksi))	GMERR(-20);
		for(k=0;k<Np;k++){
			for(j=0;j<Nst;j++){
				marray3d_set(self->Z_Ksi,t,k,j,gsl_matrix_get(self->ketatTheta_Ksi,k,j));
			}
		}
	}
	//H = Sigma0
	if(A_SET_ROW(H_Ksi,Sigma0,0))										GMERR(-22);

	//_,_,T,R,Q,a0,P0 = assemble_matrices(Np,Nstates,delta,0,sigU,sigA,sigmu,sigalph,approx)
	if(NGP_assemble_matrices(Np,Nst,self->delta,&zeroEps,self->sigKsi,self->sigA,
	&self->sigMu,&self->sigAlph,1,self->ZDummy_Ksi,self->HDummy_Ksi,
	self->T_Ksi,self->R_Ksi,self->Q_Ksi,self->a0_Ksi,self->P0_Ksi))	GMERR(-31);

	if(_NGPsample_operations(self->KsiSampler))		GMERR(-41);

	//Ksi = samp[1:Nst, :]
	if(marray3d_get_matrixSliceX(self->alphaSamples_Ksi,self->Ksi,0,0,Nst,self->work_Ksi))	GMERR(-51);

	//DKsi = samp[(Nst + 1):2Nst, :]
	if(marray3d_get_matrixSliceX(self->alphaSamples_Ksi,self->DKsi,0,Nst,2*Nst,self->work_Ksi))	GMERR(-61);
	
	//A = samp[(2Nst + 1):3Nst, :]
	if(marray3d_get_matrixSliceX(self->alphaSamples_Ksi,self->A,0,2*Nst,3*Nst,self->work_Ksi))	GMERR(-71);

	return(0);
GMERRH("NGPlaf2_sample_Ksi",1);
}

int NGPlaf2_sample_sigKsi(NGPlaf2 *self){
	double r,draw,IGa,IGb;
	double dKsit1,dKsit,at,dt,bKsi=self->bKsi;
	int i,j;
	int Nt = self->Nt;
	int Nksi = self->A->size2;
	//R2Ksi = sum((diff(DKsi,2)-A[:,1:end-1]*delta).^2,2)/delta
	//NKsi = size(A,1);
	//sigmaKsi = Array(Float64,NKsi)
	for(i=0;i<Nksi;i++){
		r = 0;
		for(j=0;j<Nt-1;j++){
            dKsit1 = gsl_matrix_get(self->DKsi,j+1,i);
            dKsit = gsl_matrix_get(self->DKsi,j,i);
            dt = gsl_vector_get(self->delta,j);
            at = gsl_matrix_get(self->A,j,i);
            r += (dKsit1-dKsit-at*dt)*(dKsit1-dKsit-at*dt)/dt;
		}
		gsl_vector_set(self->RKsi_sigKsi,i,r);
	}
	IGa = self->aKsi+(double)Nt/2;

	for(i=0;i<Nksi;i++){
		//sigKsi[i] = sqrt(rand(InverseGamma(aKsi+0.5*Nt,bKsi+0.5*R2Ksi[idx])))
        IGb = bKsi + gsl_vector_get(self->RKsi_sigKsi,i)/2.0;
        draw = sqrt(1.0/gsl_ran_gamma(self->rand,IGa,1/IGb));
        gsl_vector_set(self->sigKsi,i,draw);
	}
	
	return(0);
GMERRH("sample_sigKsi",1);
}

int NGPlaf2_sample_sigA(NGPlaf2 *self){
	double r,draw,IGa,IGb;
	double at1,at,bA;
	int i,j;
	int NA = self->A->size2;
	int Nt = self->Nt;

	//R2A = sum(diff(A,2).^2,2)/delta
	for(i=0;i<NA;i++){
		r =0;
		for(j=0;j<Nt-1;j++){
			at1 = gsl_matrix_get(self->A,j+1,i);
			at = gsl_matrix_get(self->A,j,i);
			r += (at1-at)*(at1-at)/gsl_vector_get(self->delta,j);
		}
		gsl_vector_set(self->RA_sigA,i,r);
	}
	bA = self->bA;
	IGa = self->aA + (double)Nt/2;
	for(i=0;i<NA;i++){
		IGb = bA + gsl_vector_get(self->RA_sigA,i)/2;
		//sigA[idx] = sqrt(rand(InverseGamma(aA+0.5*Nt,bA+0.5*R2A[idx])))
		draw = sqrt(1.0/gsl_ran_gamma(self->rand,IGa,1/IGb));
		gsl_vector_set(self->sigA,i,draw);
	}
	
	return(0);
GMERRH("NGPlaf2_sample_sigA",1);
}

int NGPlaf2_sample_Sigma0(NGPlaf2 *self){
	int Nt = self->Nt;
	int Np = self->Np;
	int t;
	double draw,IGa,IGb;

	for(t=0;t<Nt;t++){
		//yhat[:,t] = kron(slice(eta,:,:,t)',Theta)*slice(ksi,:,t)
		if(M_GET_ROW(eta,etat_Sigma0,t))				GMERR(-11);
		if(M_GET_ROW(Ksi,Ksit_Sigma0,t))				GMERR(-21);
		if(kron_vm(self->etat_Sigma0,self->Theta,self->kron_Sigma0))	GMERR(-31);
		if(DDOT_MV(kron_Sigma0,Ksit_Sigma0,yhatt_Sigma0))				GMERR(-41);
		if(M_SET_ROW(yhat_Sigma0,yhatt_Sigma0,t))						GMERR(-51);
	}
	//SS = sunabs2(y-yhat,2)
	if(M_SUB(yhat_Sigma0,y))							GMERR(-61);
	if(sum_abs2Col(self->yhat_Sigma0,self->SS_Sigma0))	GMERR(-71);

	IGa = self->aEps + (double)Nt/2;
	for(t=0;t<Np;t++){
		//sig2[j] = rand(InverseGamma(aEps+.5*Nt,bEps+.5*SS[j]))
		IGb = self->bEps + gsl_vector_get(self->SS_Sigma0,t)/2.0;
		draw = 1.0/gsl_ran_gamma(self->rand,IGa,1/IGb);
		gsl_matrix_set(self->Sigma0,t,t,draw);
	}
	return(0);
GMERRH("NGPlaf2_sample_Sigma0",1);
}


int NGPlaf2_sample_Psi(NGPlaf2 *self){
	int t,Nt,j,NK,Np,k,Nst,NL;
	double temp;
	double zeroEps = 0.0;
	Nt = self->Nt;
	NK = self->NK;
	Np = self->Np;
	Nst = 3*NK;
	NL = self->NL;

	//Z = zeros(Float64,Nt,Np,3*NK)
	//H = Array(Float64,Nt,Np,Np)
	marray3d_set_zero(self->Z_Psi);
	marray3d_set_zero(self->H_Psi);

	for(t=0;t<Nt;t++){
		//Z_t = Theta *reshape(slice(Ksi,:,t),NL,NK)
		if(M_GET_ROW(Ksi,Ksit_Psi,t))							GMERR(-11);
		if(reshape(self->Ksit_Psi,NL,NK,self->KsitM_Psi))		GMERR(-21);
		if(DDOT_MM(Theta,KsitM_Psi,Z_t_Psi))					GMERR(-31);

		//Z[:,1:NK,t] = Z_t
		for(k=0;k<Np;k++){
			for(j=0;j<NK;j++){
				marray3d_set(self->Z_Psi,t,k,j,gsl_matrix_get(self->Z_t_Psi,k,j));
			}
		}

		//H[:,:,t] = Z_t*Z_t'+Sigma0
		if(gTranspose(self->Z_t_Psi,self->Z_tT_Psi))			GMERR(-41);
		if(DDOT_MM(Z_t_Psi,Z_tT_Psi,ZZT_Psi))					GMERR(-51);
		if(M_ADD(ZZT_Psi,Sigma0))								GMERR(-61);
		if(marray3d_set_X(self->H_Psi,self->ZZT_Psi,t))			GMERR(-62);
	}

	if(NGP_assemble_matrices(Np,NK,self->delta,&zeroEps,self->sigPsi,self->sigB,&self->sigMu,
							&self->sigAlph,1,self->ZDummy_Psi,self->HDummy_Psi,
							self->T_Psi,self->R_Psi,self->Q_Psi,self->a0_Psi,
							self->P0_Psi))						GMERR(-71);

	//samp = NGPtools.sample(y,1,Z,H,delta,sigPsi,sigB,sigmu,sigalph,approx=true)[:,:,1]
	if(_NGPsample_operations(self->PsiSampler))					GMERR(-81);

	//Psi = samp(1:NK,:)
	if(marray3d_get_matrixSliceX(self->alphaSamples_Psi,self->Psi,0,0,NK,self->work_Psi))		GMERR(-91);

	//DPsi = samp(Nk+1:2NK,:)
	if(marray3d_get_matrixSliceX(self->alphaSamples_Psi,self->DPsi,0,NK,2*NK,self->work_Psi))	GMERR(-101);

	//B = samp[(2NK+1:3NK,:),:]
	if(marray3d_get_matrixSliceX(self->alphaSamples_Psi,self->B,0,2*NK,3*NK,self->work_Psi))	GMERR(-111);

	return(0);
GMERRH("NGPlaf2_sample_Psi",1);
}

int NGPlaf2_sample_sigPsi(NGPlaf2 *self){
	double r,draw,IGa,IGb;
	double dPsit1,dPsit,bt,dt,bPsi=self->bPsi;
	int i,j;
	int Nt = self->Nt;
	int Npsi = self->B->size2;
	for(i=0;i<Npsi;i++){
		r = 0;
		for(j=0;j<Nt-1;j++){
			dPsit1 = gsl_matrix_get(self->DPsi,j+1,i);
			dPsit = gsl_matrix_get(self->DPsi,j,i);
			dt = gsl_vector_get(self->delta,j);
			bt = gsl_matrix_get(self->B,j,i);
			r += (dPsit1-dPsit-bt*dt)*(dPsit1-dPsit-bt*dt)/dt;
		}
		gsl_vector_set(self->RPsi_sigPsi,i,r);
	}
	IGa = self->aPsi+(double)Nt/2;
	for(i=0;i<Npsi;i++){
		IGb = bPsi + gsl_vector_get(self->RPsi_sigPsi,i)/2.0;
		draw = sqrt(1.0/gsl_ran_gamma(self->rand,IGa,1/IGb));
		gsl_vector_set(self->sigPsi,i,draw);
	}

	return(0);
GMERRH("NGPlaf2_sample_sigPsi",1);
}

int NGPlaf2_sample_sigB(NGPlaf2 *self){
	double r,draw,IGa,IGb;
	double bt1,bt,bB;
	int i,j;
	int NB = self->B->size2;
	int Nt = self->Nt;

	//R2B = sum(diff(B,2).^2,2)/delta
	for(i=0;i<NB;i++){
		r =0;
		for(j=0;j<Nt-1;j++){
			bt1 = gsl_matrix_get(self->B,j+1,i);
			bt = gsl_matrix_get(self->B,j,i);
			r += (bt1-bt)*(bt1-bt)/gsl_vector_get(self->delta,j);
		}
		gsl_vector_set(self->RB_sigB,i,r);
	}

	IGa = self->aB + (double)Nt/2;
	bB = self->bB;
	for(i=0;i<NB;i++){
		IGb = bB + gsl_vector_get(self->RB_sigB,i)/2.0;
		//sigB[idx] = sqrt(rand(InverseGamma(aB+0.5*Nt,bB+0.5*R2B[idx])))
		draw = sqrt(1/gsl_ran_gamma(self->rand,IGa,1/IGb));
		gsl_vector_set(self->sigB,i,draw);
	}
	
	return(0);
GMERRH("NGPlaf2_sample_sigB",1);
}

int NGPlaf2_sample_eta(NGPlaf2 *self){
	int t,i,k;
	double temp;
	int Nt = self->Nt;
	int Np = self->Np;
	int NK = self->NK;
	int NL = self->NL;

	for(t=0;t<Nt;t++){
		//this_Ksi = reshape(slice(Ksi,:,t),NL,NK)
		if(M_GET_ROW(Ksi,Ksit_Eta,t))							GMERR(-11);
		if(M_GET_ROW(Psi,Psit_Eta,t))							GMERR(-21);
		if(reshape(self->Ksit_Eta,NL,NK,self->KsitM_Eta))		GMERR(-31);

		//Z_t = Theta *this_Ksi
		if(DDOT_MM(Theta,KsitM_Eta,Z_t_Eta))					GMERR(-41);

		//yhat = Theta * this_Ksi * slice(Psi,:,t)
		if(DDOT_MV(Z_t_Eta,Psit_Eta,yhat_Eta))					GMERR(-51);
		if(gTranspose(self->Z_t_Eta,self->Z_tT_Eta))			GMERR(-61);

		//Lambda = eye(size(Z_t,2)) + Z_t' * scale(1./diag(Sigma0),Z_t)
		for(i=0;i<Np;i++){
			for(k=0;k<NK;k++){
			//This is doing the scaling part shown above
			temp = gsl_matrix_get(self->Z_tT_Eta,k,i);
			temp = temp/gsl_matrix_get(self->Sigma0,i,i);
			gsl_matrix_set(self->Z_tT_Eta,k,i,temp);
			}
		}
		if(DDOT_MM(Z_tT_Eta,Z_t_Eta,Lambda_Eta))				GMERR(-71);
		if(M_ADD(Lambda_Eta,eye_Eta))							GMERR(-81);

		//Sigma = inv(Lambda)
		if(inv_LU(self->Lambda_Eta,self->Sigma_Eta))			GMERR(-91);

		//nu = rand(MvNormal(Sigma * scale(Z_t',1./diag(Sigma0)) * (y[:,t]-yhat),Sigma))
		//remember we scaled Z_tT above already so we didn't have to
		//make a separate temprorary variable
		if(M_GET_ROW(y,yt_Eta,t))								GMERR(-101);
		if(V_SUB(yt_Eta,yhat_Eta))								GMERR(-111);
		if(DDOT_MM(Sigma_Eta,Z_tT_Eta,SigZ_Eta))				GMERR(-121);
		if(DDOT_MV(SigZ_Eta,yt_Eta,mu_Eta))						GMERR(-131);
		if(rmvnorm(self->rand,NK,self->mu_Eta,self->Sigma_Eta,self->etat_Eta))GMERR(-141);

		//eta[:,t] = Psi[:,t] + nu
		if(V_ADD(etat_Eta,Psit_Eta))							GMERR(-151);
		if(M_SET_ROW(eta,etat_Eta,t))							GMERR(-161);
	}
	
	return(0);
GMERRH("NGPlaf2_sample_eta",1);
}

int NGPlaf2_sample_Theta(NGPlaf2 *self){
	int Nt = self->Nt;
	int Np = self->Np;
	int NL = self->NL;
	int NK = self->NK;
	int t,i;
	double sig2j;
	
	//eta_tilde = Array(Float64,NL,Nt);
	for(t=0;t<Nt;t++){
		//eta_tilde[:,t] = reshape(Ksi[:,t],NL,NK) * eta[:,t] 
		if(M_GET_ROW(Ksi,Ksit_Theta,t))										GMERR(-11);
		if(reshape(self->Ksit_Theta,NL,NK,self->KsitM_Theta))				GMERR(-21);
		if(M_GET_ROW(eta,etat_Theta,t))										GMERR(-31);
		if(DDOT_MV(KsitM_Theta,etat_Theta,eta_tildet_Theta))				GMERR(-41);
		if(M_SET_ROW(eta_tilde_Theta,eta_tildet_Theta,t))					GMERR(-51);
	}

	//etaeta = eta_tilde*eta_tilde'
	if(gTranspose(self->eta_tilde_Theta,self->eta_tildeT_Theta))			GMERR(-61);
	if(DDOT_MM(eta_tildeT_Theta,eta_tilde_Theta,etaeta_Theta))				GMERR(-71);

	for(t=0;t<Np;t++){
		//sig2j = Sigma0[,t,t]
		sig2j = gsl_matrix_get(self->Sigma0,t,t);

		//Sigmainv = (etaeta/sig2j) + diagm(squeeze(Phi[j,:],1).*tau)
		if(gsl_matrix_memcpy(self->Siginv_Theta,self->etaeta_Theta))		GMERR(-81);
		if(gsl_matrix_scale(self->Siginv_Theta,1.0/sig2j))					GMERR(-91);
		if(M_GET_ROW(Phi,Phit_Theta,t))										GMERR(-101);
		if(gsl_vector_mul(self->Phit_Theta,self->tau))						GMERR(-111);
		for(i=0;i<NL;i++){
			gsl_matrix_set(self->Phitau_Theta,i,i,gsl_vector_get(self->Phit_Theta,i));
		}
		if(gsl_matrix_add(self->Siginv_Theta,self->Phitau_Theta))			GMERR(-121);
	
		//SigmaTheta = inv(Sigmainv+tol)
		if(gsl_matrix_add(self->Siginv_Theta,self->tol_Theta))				GMERR(-131);
		if(inv_LU(self->Siginv_Theta,self->SigTheta_Theta))					GMERR(-141);

		//mu = squeeze(SigmaTheta *(eta_tilde * y[t,:]')/sig2j,2)
		if(gsl_matrix_get_col(self->yp_Theta,self->y,t))					GMERR(-151);
		if(DDOT_MV(eta_tildeT_Theta,yp_Theta,etay_Theta))					GMERR(-161);
		if(DDOT_MV(SigTheta_Theta,etay_Theta,mu_Theta))						GMERR(-171);
		if(gsl_vector_scale(self->mu_Theta,1.0/sig2j))						GMERR(-181);

		//Theta[j,:] = rand(MvNormal(mu,SigmaTheta)
		if(rmvnorm(self->rand,NL,self->mu_Theta,self->SigTheta_Theta,self->Thetat_Theta))			GMERR(-191);
		if(M_SET_ROW(Theta,Thetat_Theta,t))									GMERR(-201);
	}

	return(0);
GMERRH("NGPlaf2_sample_Theta",1);
}

int NGPlaf2_sample_Phi(NGPlaf2 *self){
	int j,l;
	int Np = self->Np;
	int NL = self->NL;
	double taul,thetajl,draw;
	for(j=0;j<Np;j++){
		for(l=0;l<NL;l++){
			taul = gsl_vector_get(self->tau,l);
			thetajl = gsl_matrix_get(self->Theta,j,l);
			thetajl = thetajl*thetajl;
			draw = gsl_ran_gamma(self->rand,2,2/(3+taul*thetajl));
			gsl_matrix_set(self->Phi,j,l,draw);
		}
	}
	return(0);
GMERRH("NGPlaf2_sample_Phi",1);
}

int NGPlaf2_sample_Tau(NGPlaf2 *self){
	int h;
	double A1,A2,a1=self->a1,a2=self->a2,draw;
	double Np = self->Np;double NL = self->NL;
	double d=0;
	double *dot=&d;

	//PhiTheta = squeeze(sum(Phi.*Theta.^2,1),1)
	gsl_matrix_memcpy(self->Theta2_Tau,self->Theta);
	gsl_matrix_mul_elements(self->Theta2_Tau,self->Theta);
	gsl_matrix_mul_elements(self->Theta2_Tau,self->Phi);
	if(sum_col(self->Theta2_Tau,self->sum_Tau))										GMERR(-11);

	//tau_minus = cumprod(theta)/theta[1]
	if(cumProd(self->theta,self->tau_minus_Tau))									GMERR(-21);
	if(gsl_vector_scale(self->tau_minus_Tau,1/gsl_vector_get(self->theta,0)))  		GMERR(-31);

	//theta[1] = rand(Gamma(a1+.5*Np,1/(1+.5*dot(tau_minus,PhiTheta))))
	A1 = a1 + (double)(Np*NL)/2.0;
	if(gsl_blas_ddot(self->tau_minus_Tau,self->sum_Tau,dot))						GMERR(-32);
	draw = gsl_ran_gamma(self->rand,A1,1/(1.0+0.5*(*dot)));
	gsl_vector_set(self->theta,0,draw);
	
	for(h=0;h<NL;h++){
		//tau_minus = cumprod(theta)/theta[h]
		if(cumProd(self->theta,self->tau_minus_Tau))								GMERR(-41);
		if(gsl_vector_scale(self->tau_minus_Tau,1/gsl_vector_get(self->theta,h)))	GMERR(-51);
		if(gsl_blas_ddot(self->tau_minus_Tau,self->sum_Tau,dot))					GMERR(-61);

		//theta[h] = rand(a2+.5*Np*(NL-h+1),1/(1+.5*dot(tau_minus,Phitheta)))
		A2 = a2 + (double)(Np*(NL-h))/2.0;
		draw = gsl_ran_gamma(self->rand,A2,1/(1.0+0.5*(*dot)));
		gsl_vector_set(self->theta,h,draw);
	}

	if(cumProd(self->theta,self->tau))												GMERR(-71);

	return(0);
GMERRH("NGPlaf2_sample_Tau",1);
}

/****************************************************
* Methods for calculating the deterministic outputs *
* However, these are not used in the operations     *
* method since the require recalculating Theta *Ksi *
* multiple times.									*
****************************************************/

int NGPlaf2_calculateMu(NGPlaf2 *self,gsl_matrix *mu){
	int Nt,Np,NK,NL,t;
	gsl_vector *ksi_this;
	gsl_vector *eta_this;
	gsl_matrix *ksiM_this;
	gsl_matrix *Thetaksi;
	gsl_vector *mut;
	//gsl_matrix *ThetaksiT;
	if((self==NULL)||(mu==NULL))				GMERR(-1);
	Nt = self->Nt;
	Np = self->Np;
	NK = self->NK;
	NL = self->NL;
	if((mu->size1!=Nt)||(mu->size2!=Np))							GMERR(-11);
	ksi_this = gsl_vector_alloc(NL*NK);
	eta_this = gsl_vector_alloc(NK);
	ksiM_this = gsl_matrix_alloc(NL,NK);
	Thetaksi = gsl_matrix_alloc(Np,NK);
	mut = gsl_vector_alloc(Np);
	//ThetaksiT = gsl_matrix_alloc(NK,Np);
	for(t=0;t<Nt;t++){
		if(gsl_matrix_get_row(ksi_this,self->Ksi,t))				GMERR(-21);
		if(gsl_matrix_get_row(eta_this,self->eta,t))				GMERR(-31);
		if(reshape(ksi_this,NL,NK,ksiM_this))						GMERR(-41);
		if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,self->Theta,ksiM_this,0.0,Thetaksi))	GMERR(-51);
		if(gsl_blas_dgemv(CblasNoTrans,1.0,Thetaksi,eta_this,0.0,mut))	GMERR(-61);
		if(gsl_matrix_set_row(mu,t,mut))							GMERR(-71);
	}

	gsl_vector_free(ksi_this);
	gsl_vector_free(eta_this);
	gsl_matrix_free(Thetaksi);
	gsl_matrix_free(ksiM_this);
	gsl_vector_free(mut);
	//gsl_matrix_free(ThetaksiT);
	return(0);
GMERRH("NGPlaf2_calculateMu",1);
}

int NGPlaf2_calculateYhat(NGPlaf2 *self,gsl_matrix *yhat){
	int Nt,Np,NK,NL,t;
	gsl_vector *ksi_this;
	gsl_vector *psi_this;
	gsl_matrix *ksiM_this;
	gsl_matrix *Thetaksi;
	gsl_vector *yhatt;
	if((self==NULL)||(yhat==NULL))				GMERR(-1);
	Nt = self->Nt;
	Np = self->Np;
	NK = self->NK;
	NL = self->NL;
	if((yhat->size1!=Nt)||(yhat->size2!=Np))	GMERR(-11);
	ksi_this = gsl_vector_alloc(NL*NK);
	psi_this = gsl_vector_alloc(NK);
	ksiM_this = gsl_matrix_alloc(NL,NK);
	Thetaksi = gsl_matrix_alloc(Np,NK);
	yhatt = gsl_vector_alloc(Np);
	for(t=0;t<Nt;t++){
		if(gsl_matrix_get_row(ksi_this,self->Ksi,t))				GMERR(-21);
		if(gsl_matrix_get_row(psi_this,self->Psi,t))				GMERR(-31);
		if(reshape(ksi_this,NL,NK,ksiM_this))						GMERR(-41);
		if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,self->Theta,ksiM_this,0,Thetaksi))	GMERR(-51);
		if(gsl_blas_dgemv(CblasNoTrans,1.0,Thetaksi,psi_this,0.0,yhatt)) GMERR(-61);
		if(gsl_matrix_set_row(yhat,t,yhatt))						GMERR(-71);
	}

	gsl_vector_free(ksi_this);
	gsl_vector_free(psi_this);
	gsl_matrix_free(ksiM_this);
	gsl_matrix_free(Thetaksi);
	gsl_vector_free(yhatt);
	return(0);
GMERRH("NGPlaf2_calculateYhat",1);
}

int NGPlaf2_calculateSigma(NGPlaf2 *self,marray3d *Sigma){
	int Nt,Np,NK,NL,t;
	gsl_vector *ksi_this;
	gsl_matrix *ksiM_this;
	gsl_matrix *Thetaksi;
	gsl_matrix *ThetaksiT;
	gsl_matrix *Sigmat;
	if((self==NULL)||(Sigma==NULL))				GMERR(-1);
	//if(marrayCheckSize(Sigma,Nt,Np,Np))			GMERR(-11);
	if(marray3d_checkSizes(Sigma,Nt,Np,Np))		GMERR(-11);
	Nt = self->Nt;
	Np = self->Np;
	NK = self->NK;
	NL = self->NL;
	ksi_this = gsl_vector_alloc(NL*NK);
	Thetaksi = gsl_matrix_alloc(Np,NK);
	ksiM_this = gsl_matrix_alloc(NL,NK);
	ThetaksiT = gsl_matrix_alloc(NK,Np);
	Sigmat = gsl_matrix_alloc(Np,Np);

	for(t=0;t<Nt;t++){
		if(gsl_matrix_get_row(ksi_this,self->Ksi,t))				GMERR(-21);
		if(reshape(ksi_this,NL,NK,ksiM_this))						GMERR(-31);
		if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,self->Theta,ksiM_this,0,Thetaksi))	GMERR(-41);
		if(gTranspose(Thetaksi,ThetaksiT))							GMERR(-51);
		if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Thetaksi,ThetaksiT,0,Sigmat))		GMERR(-61);
		if(gsl_matrix_add(Sigmat,self->Sigma0))						GMERR(-71);
		if(marray3d_set_X(Sigma,Sigmat,t))							GMERR(-81);
	//	if(marray_set_row(Sigma,Sigmat,t))							GMERR(-81);
	}

	gsl_vector_free(ksi_this);
	gsl_matrix_free(ksiM_this);
	gsl_matrix_free(Thetaksi);
	gsl_matrix_free(ThetaksiT);
	gsl_matrix_free(Sigmat);

	return(0);
GMERRH("NGPlaf2_calculateSigma",1);
}

int NGPlaf2_calculateMut(NGPlaf2 *self,int t,gsl_vector *mu){
	int Nt,Np,NK,NL;
	gsl_vector *ksi_this;
	gsl_vector *eta_this;
	gsl_matrix *ksiM_this;
	gsl_matrix *Thetaksi;
	//gsl_matrix *ThetaksiT;
	if((self==NULL)||(mu==NULL))				GMERR(-1);
	Nt = self->Nt;
	Np = self->Np;
	NK = self->NK;
	NL = self->NL;
	if((mu->size!=Np))							GMERR(-11);

	//Allocations
	ksi_this = gsl_vector_alloc(NL*NK);
	eta_this = gsl_vector_alloc(NK);
	ksiM_this = gsl_matrix_alloc(NL,NK);
	Thetaksi = gsl_matrix_alloc(Np,NK);

	//Math
	if(gsl_matrix_get_row(ksi_this,self->Ksi,t))				GMERR(-21);
	if(gsl_matrix_get_row(eta_this,self->eta,t))				GMERR(-31);
	if(reshape(ksi_this,NL,NK,ksiM_this))						GMERR(-41);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,self->Theta,ksiM_this,0.0,Thetaksi))	GMERR(-51);
	if(gsl_blas_dgemv(CblasNoTrans,1.0,Thetaksi,eta_this,0.0,mu))	GMERR(-61);

	//free
	gsl_vector_free(ksi_this);
	gsl_vector_free(eta_this);
	gsl_matrix_free(ksiM_this);
	gsl_matrix_free(Thetaksi);
 
	return(0);
GMERRH("NGPlaf2_calculateMut",1);
}

int NGPlaf2_calculateYhatt(NGPlaf2 *self,int t,gsl_vector *yhat){
	int Nt,Np,NK,NL;
	gsl_vector *ksi_this;
	gsl_vector *psi_this;
	gsl_matrix *ksiM_this;
	gsl_matrix *Thetaksi;
	if((self==NULL)||(yhat==NULL))				GMERR(-1);
	Nt = self->Nt;
	Np = self->Np;
	NK = self->NK;
	NL = self->NL;
	if((yhat->size!=Np))	GMERR(-11);

	//Allocation
	ksi_this = gsl_vector_alloc(NL*NK);
	psi_this = gsl_vector_alloc(NK);
	ksiM_this = gsl_matrix_alloc(NL,NK);
	Thetaksi = gsl_matrix_alloc(Np,NK);

	//math
	if(gsl_matrix_get_row(ksi_this,self->Ksi,t))				GMERR(-21);
	if(gsl_matrix_get_row(psi_this,self->Psi,t))				GMERR(-31);
	if(reshape(ksi_this,NL,NK,ksiM_this))						GMERR(-41);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,self->Theta,ksiM_this,0,Thetaksi))	GMERR(-51);
	if(gsl_blas_dgemv(CblasNoTrans,1.0,Thetaksi,psi_this,0.0,yhat)) GMERR(-61);

	//free variables
	gsl_vector_free(ksi_this);
	gsl_vector_free(psi_this);
	gsl_matrix_free(ksiM_this);
	gsl_matrix_free(Thetaksi);

	return(0);
GMERRH("NGPlaf2_calculateYhatt",1);
}

int NGPlaf2_calculateSigmat(NGPlaf2 *self,int t,gsl_matrix *Sigma){
	int Nt,Np,NK,NL;
	gsl_vector *ksi_this;
	gsl_matrix *ksiM_this;
	gsl_matrix *Thetaksi;
	gsl_matrix *ThetaksiT;
	gsl_matrix *Sigmat;
	if((self==NULL)||(Sigma==NULL))				GMERR(-1);
	Nt = self->Nt;
	Np = self->Np;
	NK = self->NK;
	NL = self->NL;
	if((Sigma->size1!=Np)||(Sigma->size2))		GMERR(-11);

	//allocate memory
	ksi_this = gsl_vector_alloc(NL*NK);
	Thetaksi = gsl_matrix_alloc(Np,NK);
	ksiM_this = gsl_matrix_alloc(NL,NK);
	ThetaksiT = gsl_matrix_alloc(NK,Np);
	Sigmat = gsl_matrix_alloc(Np,Np);

	//math
	if(gsl_matrix_get_row(ksi_this,self->Ksi,t))				GMERR(-21);
	if(reshape(ksi_this,NL,NK,ksiM_this))						GMERR(-31);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,self->Theta,ksiM_this,0,Thetaksi))	GMERR(-41);
	if(gTranspose(Thetaksi,ThetaksiT))							GMERR(-51);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Thetaksi,ThetaksiT,0,Sigma))		GMERR(-61);
	if(gsl_matrix_add(Sigma,self->Sigma0))						GMERR(-71);

	//free variables
	gsl_vector_free(ksi_this);
	gsl_matrix_free(ksiM_this);
	gsl_matrix_free(Thetaksi);
	gsl_matrix_free(ThetaksiT);

	return(0);
GMERRH("NGPlaf2_calculateSigmat",1);
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

/****************************************
* Class method for doing operation once *
****************************************/

int NGPlaf2_LFP(){
	
	return(0);
GMERRH("NGPlaf2_LFP",1);
}

/*********************************
* Writing method outputs to file *
*********************************/

//outputs
int NGPlaf2_writeKsiOut(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	//if(marray_fprintf(fp,self->Ksi_out,"%.15f"))	GMERR(-11);
	if(marray3d_write(fp,self->Ksi_out))GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writeKsiOut",1);
}

int NGPlaf2_writePsiOut(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(marray3d_write(fp,self->Psi_out))GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writePsiOut",1);
}

int NGPlaf2_writeyhatOut(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(marray3d_write(fp,self->yhat_out))GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writeyhatOut",1);
}

int NGPlaf2_writeMuOut(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(marray3d_write(fp,self->mu_out))	GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writeMuOut",1);
}

int NGPlaf2_writeSigmaOut(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(marray4d_write(fp,self->Sigma_out))GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writeSigmaOut",1);
}

int NGPlaf2_writeOutputs(NGPlaf2 *self,char *baseName){
	//Here we have to make all the output names
	FILE *fp;
	char strK[16] = "_Ksi.txt";
	char strP[16] = "_Psi.txt";
	char strY[16] = "_yhat.txt";
	char strM[16] = "_Mu.txt";
	char strS[16] = "_Sigma0.txt";
	char outK[1024],outP[1024],outY[1024],outM[1024],outS[1024];
	sprintf(outK,"%s%s",baseName,strK);
	sprintf(outP,"%s%s",baseName,strP);
	sprintf(outY,"%s%s",baseName,strY);
	sprintf(outM,"%s%s",baseName,strM);
	sprintf(outS,"%s%s",baseName,strS);

	//Here we actually write out all the data
	if(self==NULL)							GMERR(-1);
	if(baseName==NULL)						GMERR(-3);
	if(NGPlaf2_writeKsiOut(self,outK))		GMERR(-11);
	if(NGPlaf2_writePsiOut(self,outP))		GMERR(-21);
	if(NGPlaf2_writeyhatOut(self,outY))		GMERR(-31);
	if(NGPlaf2_writeMuOut(self,outM))		GMERR(-41);
	if(NGPlaf2_writeSigmaOut(self,outS))		GMERR(-51);
	
	return(0);
GMERRH("NGPlaf2_writeOutputs",1);
}


//Individual samples
int NGPlaf2_writeKsi(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_matrix_fprintf(fp,self->Ksi,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writeKsi",1);
}

int NGPlaf2_writesigKsi(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_vector_fprintf(fp,self->sigKsi,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writesigKsi",1);
}

int NGPlaf2_writesigA(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_vector_fprintf(fp,self->sigA,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writesigA",1);
}

int NGPlaf2_writeSigma0(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_matrix_fprintf(fp,self->Sigma0,"%.15f"))	GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writeSigma0",1);
}

int NGPlaf2_writePsi(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_matrix_fprintf(fp,self->Psi,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writePsi",1);
}

int NGPlaf2_writesigPsi(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_vector_fprintf(fp,self->sigPsi,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writesigPsi",1);
}

int NGPlaf2_writesigB(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_vector_fprintf(fp,self->sigB,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writesigB",1);
}

int NGPlaf2_writeeta(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_matrix_fprintf(fp,self->eta,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writeeta",1);
}

int NGPlaf2_writeTheta(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_matrix_fprintf(fp,self->Theta,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writeTheta",1);
}

int NGPlaf2_writePhi(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_matrix_fprintf(fp,self->Phi,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writePhi",1);
}

int NGPlaf2_writetau(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_vector_fprintf(fp,self->tau,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writeSigmaOut",1);
}

int NGPlaf2_writetheta(NGPlaf2 *self,char *name){
	FILE *fp;
	if(self==NULL)						GMERR(-1);
	fp = fopen(name,"w");
	if(fp==NULL)						GMERR(-3);
	if(gsl_vector_fprintf(fp,self->theta,"%.15f"))		GMERR(-11);
	fclose(fp);
	return(0);
GMERRH("NGPlaf2_writetheta",1);
}

/**********************************************
* Here we have the methods for printing some  *
* of the variables in the object. This won't  *
* really be used once the code is correct but *
* its helpful when debugging and I might as	  *
* well leave them in.						  *
**********************************************/
int NGPlaf2_printKsi(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("Ksi\n");
	printGSLMatrix(self->Ksi);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printKsi",1);
}

int NGPlaf2_printsigKsi(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("sigKsi\n");
	printGSLVectorT(self->sigKsi);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printsigKsi",1);
}

int NGPlaf2_printsigA(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("A\n");
	printGSLVectorT(self->sigA);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printsigA",1);
}

int NGPlaf2_printSigma0(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("Sigma0\n");
	printGSLMatrix(self->Sigma0);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printSigma0",1);
}

int NGPlaf2_printPsi(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("Psi\n");
	printGSLMatrix(self->Psi);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printPsi",1);
}

int NGPlaf2_printsigPsi(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("sigPsi\n");
	printGSLVectorT(self->sigPsi);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printsigPsi",1);
}

int NGPlaf2_printsigB(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("B\n");
	printGSLVectorT(self->sigB);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printsigB",1);
}

int NGPlaf2_printtheta(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("theta\n");
	printGSLVectorT(self->theta);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printtheta",1);
}

int NGPlaf2_printTheta(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("Theta\n");
	printGSLMatrix(self->Theta);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printTheta",1);
}

int NGPlaf2_printPhi(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("Phi\n");
	printGSLMatrix(self->Phi);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printPhi",1);
}

int NGPlaf2_printtau(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("tau\n");
	printGSLVectorT(self->tau);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printtau",1);
}

int NGPlaf2_printeta(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("eta\n");
	printGSLMatrix(self->eta);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printeta",1);
}

int NGPlaf2_printDKsi(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("DKsi\n");
	printGSLMatrix(self->DKsi);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printDKsi",1);
}

int NGPlaf2_printA(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("A\n");
	printGSLMatrix(self->A);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printA",1);
}

int NGPlaf2_printDPsi(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("DPsi\n");
	printGSLMatrix(self->DPsi);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printDPsi",1);
}

int NGPlaf2_printB(NGPlaf2 *self){
	if(self==NULL)				GMERR(-1);
	printf("\n");
	printf("B\n");
	printGSLMatrix(self->B);
	printf("\n");
	return(0);
GMERRH("NGPlaf2_printB",1);
}

