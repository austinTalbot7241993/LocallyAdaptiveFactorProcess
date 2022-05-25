#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include "gmlib.h"
#include "matrix_utils.h"
#include "mdarray.h"
#include "NGPmcmc.h"

int main(){
	FILE *fp;
//	int nSamps = 10000;
	int i;
	int Niter = 2000;
	int Nt = 1001;
//	gsl_vector *epsSamps = gsl_vector_alloc(nSamps);
	gsl_matrix *y = gsl_matrix_alloc(Nt,1);
	gsl_vector *sigU = gsl_vector_alloc(1);
	gsl_vector *sigA = gsl_vector_alloc(1);
	double a = 1.0;
	double b = 1.0;
	double sigMu = 4.0;
	double sigAlph = 4.0;
	double sigEps = 2.0;
	marray3d *th = marray3d_alloc(Niter,Nt,3);
	gsl_matrix *sig = gsl_matrix_alloc(Niter,3);
	gsl_vector *tobs = gsl_vector_alloc(Nt);
	NGPmcmc *myMCMC;

	gsl_vector_set(sigU,0,1000.0);
	gsl_vector_set(sigA,0,5);

	printf("ahd;lfkjda;fl\n");
	fp = fopen("y.txt","r");
	if(fp==NULL)			GMERR(-11);
	gsl_matrix_fscanf(fp,y);
	fclose(fp);
	printf("ahd;lfkjda;fl\n");

	fp = fopen("tobs.txt","r");
	if(fp==NULL)			GMERR(-21);
	gsl_vector_fscanf(fp,tobs);
	fclose(fp);

	printf("Hello\n");
	myMCMC = NGPmcmc_New();
	printf("Hello2\n");
	if(NGPmcmc_construct(myMCMC,y,tobs,Niter,sigU,sigA,sigEps,sigMu,sigAlph,a,b,th,sig))	GMERR(-31);

	if(NGPmcmc_operations(myMCMC))			GMERR(-41);
	if(NGPmcmc_writeOutputs(myMCMC,"AMY_M1"))	GMERR(-101);

//	fp = fopen("U.txt","r");
//	if(fp==NULL)		GMERR(-31);
//	gsl_vector_fscanf(fp,myMCMC->U);
//	printGSLVector(myMCMC->U);
//	printf(">>>>>>\n");
//	fclose(fp);
//	for(i=0;i<nSamps;i++){
//		if(NGPmcmc_drawSigEps(myMCMC))		GMERR(-41);
//		gsl_vector_set(epsSamps,i,myMCMC->sigEps);
//	}
//	fp = fopen("epsSamps.txt","w");
//	if(fp==NULL)			GMERR(-51);
//	gsl_vector_fprintf(fp,epsSamps,"%0.15f");
//	fclose(fp);


//	fp = fopen("G.txt","r");
//	if(fp==NULL)				GMERR(-31);
//	if(marray3d_read(fp,myMCMC->G))	GMERR(-32);
//	fclose(fp);
//
//	fp = fopen("Gtilde.txt","r");
//	if(fp==NULL)				GMERR(-41);
//	if(marray3d_read(fp,myMCMC->Gtilde))	GMERR(-42);
//	fclose(fp);
//
//	fp = fopen("W.txt","r");
//	if(fp==NULL)				GMERR(-61);
//	if(marray3d_read(fp,myMCMC->W))	GMERR(-62);
//	fclose(fp);
//
//	fp = fopen("Htilde.txt","r");
//	if(fp==NULL)				GMERR(-71);
//	if(marray3d_read(fp,myMCMC->Htilde))	GMERR(-72);
//	fclose(fp);
//
//	fp = fopen("Wtilde.txt","r");
//	if(fp==NULL)				GMERR(-81);
//	if(marray3d_read(fp,myMCMC->Wtilde))	GMERR(-82);
//	fclose(fp);
//
//	fp = fopen("Wstar.txt","r");
//	if(fp==NULL)				GMERR(-91);
//	if(marray3d_read(fp,myMCMC->Wstar))	GMERR(-92);
//	fclose(fp);
//	
//	fp = fopen("Wtilde_star.txt","r");
//	if(fp==NULL)				GMERR(-91);
//	if(marray3d_read(fp,myMCMC->Wtildestar))	GMERR(-92);
//	fclose(fp);
//
//	fp = fopen("theta.txt","r");
//	if(fp==NULL)				GMERR(-101);
//	gsl_matrix_fscanf(fp,myMCMC->theta);
//	fclose(fp);
//	
//	fp = fopen("theta_star.txt","r");
//	if(fp==NULL)				GMERR(-111);
//	gsl_matrix_fscanf(fp,myMCMC->thetastar);
//	fclose(fp);
//
//	if(NGPmcmc_calc_acceptance_prob(myMCMC))	GMERR(-221);

//	if(NGPmcmc_writeOutputs(myMCMC,"test2"))	GMERR(-101);
	

	return(0);
GMERRH("main",1);
}
