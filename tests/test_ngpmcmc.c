#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <lafp/gmlib.h>
#include <lafp/matrix_utils.h>
#include <lafp/mdarray.h>
#include <lafp/NGPmcmc.h>

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

	NGPmcmc_free(myMCMC);
	gsl_matrix_free(y);
	gsl_vector_free(sigU);
	gsl_vector_free(sigA);
	marray3d_free(th);
	gsl_matrix_free(sig);
	gsl_vector_free(tobs);

	return(0);
GMERRH("main",1);
}
