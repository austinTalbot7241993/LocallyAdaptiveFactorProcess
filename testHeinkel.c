#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include "gmlib.h"
#include "matrix_utils.h"

int main(){
	int i,N,N2;
	double temp;
	N = 3;
	N2 = 5;
	gsl_vector *c = gsl_vector_alloc(N);
	gsl_vector *r = gsl_vector_alloc(N2);
	gsl_matrix *out = gsl_matrix_alloc(N,N2);
	for(i=0;i<N;i++){
		temp = 2.1*i;
		gsl_vector_set(c,i,temp);
	}
	gsl_vector_set(r,0,10.0);
	gsl_vector_set(r,1,11.0);
	gsl_vector_set(r,2,12.0);
	gsl_vector_set(r,3,14.0);
	gsl_vector_set(r,4,15.0);

	if(hankel(c,r,out))		GMERR(-11);
	printGSLMatrix(out);
	return(0);
GMERRH("main",1);
}
