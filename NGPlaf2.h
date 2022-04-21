#ifndef _ABT_NGPLAF2_H
#define _ABT_NGPLAF2_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include "NGPsample2.h"
#include "mdarray.h"

typedef struct{
	int Np;
	int Nt;
	int NK;
	int NL;
	int Nst;//Number of sttates in dynamical model
	int Nm;//Number of states in nGP exanded model

	//Priors
	double sigMu;
	double sigAlph;
	double aEps,bEps;
	double aKsi,bKsi;
	double aA,bA;
	double aPsi,bPsi;
	double aB,bB;
	double a1,a2;

	int Niter,burnin,Nthin;

	//Inputs
	gsl_matrix *y;
	gsl_vector *tobs;
	gsl_vector *delta;

	/**********
	* outputs *
	**********/

	marray3d *Ksi_out;
	marray3d *Psi_out;
	marray3d *yhat_out;
	marray3d *mu_out;
	marray4d *Sigma_out;

	/********************
	* overall variables *
	********************/
	gsl_matrix *eta;
	gsl_matrix *Theta;
	gsl_matrix *Phi;
	gsl_vector *theta;
	gsl_vector *tau;
	gsl_matrix *Sigma0;
	gsl_matrix *Ksi;
	gsl_vector *sigKsi;
	gsl_vector *sigA;
	gsl_matrix *Psi;
	gsl_vector *sigPsi;
	gsl_vector *sigB;
	gsl_matrix *DKsi;
	gsl_matrix *DPsi;
	gsl_matrix *A;
	gsl_matrix *B;

	/***************
	* sample_Theta *
	***************/
	gsl_vector *eta_tildet_Theta;
	gsl_vector *etat_Theta;
	gsl_vector *etay_Theta;
	gsl_vector *Ksit_Theta;
	gsl_vector *mu_Theta;
	gsl_vector *Phit_Theta;
	gsl_vector *Thetat_Theta;
	gsl_vector *yp_Theta;

	gsl_matrix *eta_tilde_Theta;
	gsl_matrix *eta_tildeT_Theta;
	gsl_matrix *etaeta_Theta;
	gsl_matrix *KsitM_Theta;
	gsl_matrix *Phitau_Theta;
	gsl_matrix *Siginv_Theta;
	gsl_matrix *SigTheta_Theta;
	gsl_matrix *tol_Theta;

	/*************
	* sample_Ksi *
	*************/
	marray3d *ZDummy_Ksi;
	marray3d *HDummy_Ksi;
	gsl_vector *etat_Ksi;
	gsl_matrix *ketatTheta_Ksi;

	//all used for the sampler
	marray3d *Z_Ksi;
	marray3d *H_Ksi;
	marray3d *T_Ksi;
	marray3d *R_Ksi;
	marray3d *Q_Ksi;
	marray3d *alphaSamples_Ksi;
	gsl_matrix *P0_Ksi;
	gsl_vector *a0_Ksi;
	gsl_vector *work_Ksi;
	_NGPsample *KsiSampler;

	/****************
	* sample_sigKsi *
	****************/
	gsl_vector *RKsi_sigKsi;

	/**************
	* sample_sigA *
	**************/
	gsl_vector *RA_sigA;

	/*************
	* sample_Psi *
	*************/
	marray3d *ZDummy_Psi;
	marray3d *HDummy_Psi;
	gsl_vector *Ksit_Psi;
	gsl_matrix *KsitM_Psi;
	gsl_matrix *Z_t_Psi;
	gsl_matrix *Z_tT_Psi;
	gsl_matrix *ZZT_Psi;

	//all used for the sampler
	marray3d *Z_Psi;
	marray3d *H_Psi;
	marray3d *T_Psi;
	marray3d *R_Psi;
	marray3d *Q_Psi;
	marray3d *alphaSamples_Psi;
	gsl_matrix *P0_Psi;
	gsl_vector *a0_Psi;
	gsl_vector *work_Psi;
	_NGPsample *PsiSampler;

	/****************
	* sample_sigPsi *
	****************/
	gsl_vector *RPsi_sigPsi;

	/**************
	* sample_sigB *
	**************/
	gsl_vector *RB_sigB;

	/*************
	* sample_eta *
	*************/
	gsl_vector *etat_Eta;
	gsl_vector *Ksit_Eta;
	gsl_vector *mu_Eta;
	gsl_vector *Psit_Eta;
	gsl_vector *yhat_Eta;
	gsl_vector *yt_Eta;

	gsl_matrix *eye_Eta;
	gsl_matrix *KsitM_Eta;
	gsl_matrix *Lambda_Eta;
	gsl_matrix *Sigma_Eta;
	gsl_matrix *SigZ_Eta;
	gsl_matrix *Z_t_Eta;
	gsl_matrix *Z_tT_Eta;

	/****************
	* sample_Sigma0 *
	****************/
	gsl_vector *etat_Sigma0;
	gsl_vector *Ksit_Sigma0;
	gsl_vector *yhatt_Sigma0;
	gsl_vector *SS_Sigma0;

	gsl_matrix *kron_Sigma0;
	gsl_matrix *yhat_Sigma0;


	/*************
	* sample_tau *
	*************/
	gsl_matrix *Theta2_Tau;
	gsl_vector *sum_Tau;
	gsl_vector *tau_minus_Tau;

	/*************
	* operations *
	*************/
	gsl_vector *eta_t_O;
	gsl_vector *Ksi_t_O;
	gsl_vector *mu_t_O;
	gsl_vector *Psi_t_O;
	gsl_vector *yhat_t_O;

	gsl_matrix *Ksi_this_O;
	gsl_matrix *ThetaKsi_O;
	gsl_matrix *ThetaKsiT_O;
	gsl_matrix *TKKT_O;

	gsl_rng *rand;

} NGPlaf2;


/***************************************
* Constructor, destructor, main method *
***************************************/
NGPlaf2 * NGPlaf2_New();
int NGPlaf2_init(NGPlaf2 *self);
int NGPlaf2_free(NGPlaf2 *s);
int NGPlaf2_construct(NGPlaf2 *self,gsl_vector *tobs,gsl_matrix *y,int NK,int NL,int Niter,marray3d *Ksi_out,marray3d *Psi_out,marray3d *yhat_out,marray3d *mu_out,marray4d *Sigma_out,gsl_vector *sigPrior,gsl_vector *epsPrior,gsl_vector *ksiPrior,gsl_vector *APrior,gsl_vector *psiPrior,gsl_vector *BPrior,gsl_vector *aaPrior);
int NGPlaf2_operations(NGPlaf2 *self);

/****************************************
* Class method for doing operation once *
****************************************/

int NGPlaf2_LFP();

/************************************************
* Modifying priors or reset the outputs to zero *
************************************************/

int NGPlaf2_modifyPrior(NGPlaf2 *self,char parameter,gsl_vector *param);
int NGPlaf2_reset(NGPlaf2 *self);

/******************************
* Individual sampling methods *
******************************/

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

/*Calculate deterministic outputs */
int NGPlaf2_calculateMu(NGPlaf2 *self,gsl_matrix *mu);
int NGPlaf2_calculateYhat(NGPlaf2 *self,gsl_matrix *yhat);
int NGPlaf2_calculateSigma(NGPlaf2 *self,marray3d *Sigma);

int NGPlaf2_calculateMut(NGPlaf2 *self,int t,gsl_vector *mu);
int NGPlaf2_calculateYhatt(NGPlaf2 *self,int t,gsl_vector *yhat);
int NGPlaf2_calculateSigmat(NGPlaf2 *self,int t,gsl_matrix *Sigma);


/*********************************
* Writing method outputs to file *
*********************************/

//outputs
int NGPlaf2_writeKsiOut(NGPlaf2 *self,char *name);
int NGPlaf2_writePsiOut(NGPlaf2 *self,char *name);
int NGPlaf2_writeyhatOut(NGPlaf2 *self,char *name);
int NGPlaf2_writeMuOut(NGPlaf2 *self,char *name);
int NGPlaf2_writeSigmaOut(NGPlaf2 *self,char *name);
int NGPlaf2_writeOutputs(NGPlaf2 *self,char *baseName);

//Individual samples
int NGPlaf2_writeKsi(NGPlaf2 *self,char *name);
int NGPlaf2_writesigKsi(NGPlaf2 *self,char *name);
int NGPlaf2_writesigA(NGPlaf2 *self,char *name);
int NGPlaf2_writeSigma0(NGPlaf2 *self,char *name);
int NGPlaf2_writePsi(NGPlaf2 *self,char *name);
int NGPlaf2_writesigPsi(NGPlaf2 *self,char *name);
int NGPlaf2_writesigB(NGPlaf2 *self,char *name);
int NGPlaf2_writeeta(NGPlaf2 *self,char *name);
int NGPlaf2_writeTheta(NGPlaf2 *self,char *name);
int NGPlaf2_writePhi(NGPlaf2 *self,char *name);
int NGPlaf2_writetau(NGPlaf2 *self,char *name);
int NGPlaf2_writetheta(NGPlaf2 *self,char *name);

/*******************************
* Printing methods (debugging) *
*******************************/
int NGPlaf2_printKsi(NGPlaf2 *self);
int NGPlaf2_printsigKsi(NGPlaf2 *self);
int NGPlaf2_printsigA(NGPlaf2 *self);
int NGPlaf2_printSigma0(NGPlaf2 *self);
int NGPlaf2_printPsi(NGPlaf2 *self);
int NGPlaf2_printsigPsi(NGPlaf2 *self);
int NGPlaf2_printsigB(NGPlaf2 *self);
int NGPlaf2_printtheta(NGPlaf2 *self);
int NGPlaf2_printTheta(NGPlaf2 *self);
int NGPlaf2_printPhi(NGPlaf2 *self);
int NGPlaf2_printtau(NGPlaf2 *self);
int NGPlaf2_printeta(NGPlaf2 *self);
int NGPlaf2_printDKsi(NGPlaf2 *self);
int NGPlaf2_printA(NGPlaf2 *self);
int NGPlaf2_printDPsi(NGPlaf2 *self);
int NGPlaf2_printB(NGPlaf2 *self);


#endif
