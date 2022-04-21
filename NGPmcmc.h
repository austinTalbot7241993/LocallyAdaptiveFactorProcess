#ifndef _ABT_NGPMCMC_H
#define _ABT_NGPMCMC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include "gmlib.h"
#include "mdarray.h"
#include "SSsimulate2.h"

/*Code by Austin Bryan Talbot
 *        a.talbot@verizon.net
 *
 *Description:
 *
 *Completed:
 *
 *Revision History:
 *
*/

typedef struct{

    /*******************
    ********************
    ** method inputs  **
    ********************
    *******************/
	gsl_matrix *y;
	gsl_vector *tobs;
	gsl_vector *delta;

	double sigEps;
	gsl_vector *sigU;
	gsl_vector *sigA;

    /*******************
    ********************
    ** method outputs **
    ********************
    *******************/
	int Niter;
	int burn;
	int Nt;
	marray3d *th;
	gsl_matrix *sig;

	/*********
	* priors *
	*********/
	double a;
	double b;
	double sigMu;
	double sigAlph;

    /*******************************
    ********************************
    ** method temporary variables **
    ********************************
    *******************************/
	gsl_vector *y_vec;

	/*************
	* draw_state *
	*************/
	marray3d *Z_s;
	marray3d *H_s;
	marray3d *T_s;
	marray3d *R_s;
	marray3d *Q_s;
	gsl_vector *a0_s;
	gsl_matrix *P0_s;
	SSsimulate2 *drawState;

	/**************
	* draw_sigEps *
	**************/
	gsl_vector *U;

	/***********************
	* calc_acceptance_prob *
	***********************/
	double prob;
	gsl_vector *v_t;
	gsl_vector *vstar_t;
	gsl_vector *v_03;
	gsl_vector *v_02;
	gsl_vector *work3;
	gsl_vector *work2;
	gsl_vector *theta_t;
	gsl_vector *thetastar_t;
	gsl_vector *Gthetastar;

	gsl_matrix *HtildeT;

	/************
	* draw_sigs *
	************/

	gsl_vector *sigUstar;
	gsl_vector *sigAstar;

	gsl_matrix *theta;
	gsl_matrix *thetastar;

	gsl_vector *Ustar;
	gsl_vector *DUstar;
	gsl_vector *Astar;

	//Assembling the matrices

	//Regular
	marray3d *Zdummy;
	marray3d *Hdummy;
	marray3d *G;
	marray3d *H;
	marray3d *W;
	gsl_vector *a0dummy;
	gsl_matrix *P0dummy;

	//Tilde
	marray3d *Ztdummy;
	marray3d *Htdummy;
	marray3d *Gtilde;
	marray3d *Htilde;
	marray3d *Wtilde;
	gsl_vector *a0tildedummy;
	gsl_matrix *P0tildedummy;

	//Star
	marray3d *Zstardummy;
	marray3d *Hstardummy;
	marray3d *Gstardummy;
	marray3d *Tstardummy;
	marray3d *Wstar;
	gsl_vector *a0stardummy;
	gsl_matrix *P0stardummy;

	//Tilde star
	marray3d *Ztildestardummy;
	marray3d *Htildestardummy;
	marray3d *Gtildestardummy;
	marray3d *Ttildestardummy;
	marray3d *Wtildestar;
	gsl_vector *a0tildestardummy;
	gsl_matrix *P0tildestardummy;

	gsl_vector *a0_star;
	gsl_matrix *P0_star;
	marray3d *Z_star;
	marray3d *H_star;
	marray3d *T_star;
	marray3d *R_star;
	marray3d *Q_star;

	SSsimulate2 *drawThetaStar;


	/**************************
	* random number generator *
	**************************/
	gsl_rng *rand;

} NGPmcmc;

/***************************************
* Constructor, destructor, main method *
***************************************/

NGPmcmc * NGPmcmc_New();
int NGPmcmc_init(NGPmcmc *self);
int NGPmcmc_free(NGPmcmc *s);
int NGPmcmc_construct(NGPmcmc *self,gsl_matrix *y,gsl_vector *tobs,int Niter,gsl_vector *sigU,gsl_vector *sigA,double sigEps,double sigMu,double sigAlph,double a,double b,marray3d *th,gsl_matrix *sig);
int NGPmcmc_operations(NGPmcmc *self);

/*************************************
* Sampling methods within operations *
*************************************/

int NGPmcmc_drawState(NGPmcmc *self);
int NGPmcmc_drawSigEps(NGPmcmc *self);
int NGPmcmc_calc_acceptance_prob(NGPmcmc *self);
int NGPmcmc_drawSigs(NGPmcmc *self);

/****************************************
* Class method for doing operation once *
****************************************/

int NGPmcmc_NGPmcmc();

/*********************************
* Parameter modification methods *
*********************************/

int NGPmcmc_updatePrior(NGPmcmc *self,char parameter,gsl_vector *param);

/*********************************
* Writing method outputs to file *
*********************************/

int NGPmcmc_writeTheta(NGPmcmc *self,char *name);
int NGPmcmc_writeSig(NGPmcmc *self,char *name);
int NGPmcmc_writeOutputs(NGPmcmc *self,char *baseName);


#endif
