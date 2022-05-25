/***************************************************************************************
 *  Multivariate Normal density function and random number generator
 *  Multivariate Student t density function and random number generator
 *  Wishart random number generator
 *  Using GSL -> www.gnu.org/software/gsl
 *
 *  Copyright (C) 2006  Ralph dos Santos Silva
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  AUTHOR 
 *     Ralph dos Santos Silva,  address@hidden
 *     March, 2006       
***************************************************************************************/
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include "gmlib.h"
#include "mvnorm.h"
#include "matrix_utils.h"

int _gm_errNum;

/*****************************************************************************************************************/
/*****************************************************************************************************************/
int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result){
/* multivariate normal distribution random number generator */
/*
*	n	dimension of the random vetor
*	mean	vector of means of size n
*	var	variance matrix of dimension n x n
*	result	output variable with a sigle random vector normal distribution generation
*/
	int k;
	gsl_matrix *work = gsl_matrix_alloc(n,n);

	gsl_matrix_memcpy(work,var);
	if(cmatrixIsZero(var))								GMERR(-1);	
	if(positiveDefinite(var)){
		gsl_vector *eval = gsl_vector_alloc(n);
		gsl_matrix *mcpy = gsl_matrix_alloc(n,n);
		gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(n);
		gsl_matrix_memcpy(mcpy,var);
		gsl_eigen_symm(mcpy,eval,w);
		printf("\n\n");
		printGSLVectorT(eval);
		printf("\n\n");
		printGSLMatrix(var);
		gsl_vector_free(eval);
		gsl_matrix_free(mcpy);
		gsl_eigen_symm_free(w);
		GMERR(-11);
	}
	if(gsl_linalg_cholesky_decomp(work))				GMERR(-11);

	for(k=0; k<n; k++)
		gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
	gsl_vector_add(result,mean);

	gsl_matrix_free(work);

	return(0);
GMERRH("rmvnorm",1);
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var){
/* multivariate normal density function    */
/*
*	n	dimension of the random vetor
*	mean	vector of means of size n
*	var	variance matrix of dimension n x n
*/
int s;
double ax,ay;
gsl_vector *ym, *xm;
gsl_matrix *work = gsl_matrix_alloc(n,n), 
           *winv = gsl_matrix_alloc(n,n);
gsl_permutation *p = gsl_permutation_alloc(n);

gsl_matrix_memcpy( work, var );
gsl_linalg_LU_decomp( work, p, &s );
gsl_linalg_LU_invert( work, p, winv );
ax = gsl_linalg_LU_det( work, s );
gsl_matrix_free( work );
gsl_permutation_free( p );

xm = gsl_vector_alloc(n);
gsl_vector_memcpy( xm, x);
gsl_vector_sub( xm, mean );
ym = gsl_vector_alloc(n);
gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
gsl_matrix_free( winv );
gsl_blas_ddot( xm, ym, &ay);
gsl_vector_free(xm);
gsl_vector_free(ym);
ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );

return ay;
}
int ldmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var,gsl_matrix *work,gsl_matrix *winv,double *ans){
	
	int s;
	double ax,ay;
	gsl_vector *ym, *xm;
	gsl_permutation *p = gsl_permutation_alloc(n);
	gsl_matrix_memcpy( work, var );
	gsl_linalg_LU_decomp( work, p, &s );
	gsl_linalg_LU_invert( work, p, winv );
	ax = gsl_linalg_LU_det( work, s );
	gsl_matrix_free( work );
	gsl_permutation_free( p );
	xm = gsl_vector_alloc(n);
	gsl_vector_memcpy( xm, x);
	gsl_vector_sub( xm, mean );
	ym = gsl_vector_alloc(n);
	gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
	gsl_matrix_free( winv );
	gsl_blas_ddot( xm, ym, &ay);
	gsl_vector_free(xm);
	gsl_vector_free(ym);
//	ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );
	ay = (double)n*-0.5*log(2*M_PI)-0.5*log(fabs(ay));
	printf("%f\n",ay);
	ay = ay - n*0.5*log(ax);
	*ans = ay;
	printf("%f\n",ay);
	printf("%f\n",ax);
	return(0);
GMERRH("ldmvnorm",1);
}
/*
int ldmvnorm(int n,const gsl_vector *x,const gsl_vector *mean,const gsl_matrix *var,gsl_matrix *work,gsl_matrix *winv,double *ans){
	int s;
	double det,density = 0,prod=0;
	gsl_vector *xm = gsl_vector_alloc(n);
	gsl_vector *ym = gsl_vector_alloc(n);
	gsl_permutation *p = gsl_permutation_alloc(n);

	gsl_matrix_memcpy(work,var);
	gsl_linalg_LU_decomp(work,p,&s);
	gsl_linalg_LU_invert(work,p,winv);
	det = gsl_linalg_LU_det(work,s);
	gsl_permutation_free(p);

	gsl_vector_memcpy(xm,x);
	gsl_vector_sub(xm,mean);
	printf(">>>>>>>>>>>>1\n");
	printGSLVectorT(xm);
	printf(">>>>>>>>>>>>1\n");
	
	gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
	printGSLVectorT(ym);
	printGSLMatrix(winv);
	printf("\n");
	gsl_blas_ddot(xm,ym,&prod);


	density = (double)n/2.0*log(2*M_PI);
	printf("%f\n",density);
	density = -1.0*density;
	printf("%f\n",density);
	density = density - 1.0/2.0*log(fabs(det));
	printf("%f\n",density);
	density = density-0.5*prod;
	printf("%f\n",density);
	*ans = density;

	return(0);
GMERRH("ldmvnorm",1);
}
*/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
int rmvt(const gsl_rng *r, const int n, const gsl_vector *location, const gsl_matrix *scale, const int dof, gsl_vector *result){
/* multivariate Student t distribution random number generator */
/*
*	n	 dimension of the random vetor
*	location vector of locations of size n
*	scale	 scale matrix of dimension n x n
*	dof	 degrees of freedom
*	result	 output variable with a single random vector normal distribution generation
*/
int k;
gsl_matrix *work = gsl_matrix_alloc(n,n);
double ax = 0.5*dof; 

ax = gsl_ran_gamma(r,ax,(1/ax));     /* gamma distribution */

gsl_matrix_memcpy(work,scale);
gsl_matrix_scale(work,(1/ax));       /* scaling the matrix */
gsl_linalg_cholesky_decomp(work);

for(k=0; k<n; k++)
	gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
gsl_vector_add(result, location);

gsl_matrix_free(work);

return 0;
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
double dmvt(const int n, const gsl_vector *x, const gsl_vector *location, const gsl_matrix *scale, const int dof){
/* multivariate Student t density function */
/*
*	n	 dimension of the random vetor
*	location vector of locations of size n
*	scale	 scale matrix of dimension n x n
*	dof	 degrees of freedom
*/
int s;
double ax,ay,az=0.5*(dof + n);
gsl_vector *ym, *xm;
gsl_matrix *work = gsl_matrix_alloc(n,n), 
           *winv = gsl_matrix_alloc(n,n);
gsl_permutation *p = gsl_permutation_alloc(n);

gsl_matrix_memcpy( work, scale );
gsl_linalg_LU_decomp( work, p, &s );
gsl_linalg_LU_invert( work, p, winv );
ax = gsl_linalg_LU_det( work, s );
gsl_matrix_free( work );
gsl_permutation_free( p );

xm = gsl_vector_alloc(n);
gsl_vector_memcpy( xm, x);
gsl_vector_sub( xm, location );
ym = gsl_vector_alloc(n);
gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
gsl_matrix_free( winv );
gsl_blas_ddot( xm, ym, &ay);
gsl_vector_free(xm);
gsl_vector_free(ym);

ay = pow((1+ay/dof),-az)*gsl_sf_gamma(az)/(gsl_sf_gamma(0.5*dof)*sqrt( pow((dof*M_PI),n)*ax ));

return ay;
}
/*****************************************************************************************************************/
/*****************************************************************************************************************/
int rwishart(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, gsl_matrix *result){
/* Wishart distribution random number generator */
/*
*	n	 gives the dimension of the random matrix
*	dof	 degrees of freedom
*	scale	 scale matrix of dimension n x n
*	result	 output variable with a single random matrix Wishart distribution generation
*/
int k,l;
gsl_matrix *work = gsl_matrix_calloc(n,n);

for(k=0; k<n; k++){
	gsl_matrix_set( work, k, k, sqrt( gsl_ran_chisq( r, (dof-k) ) ) );
	for(l=0; l<k; l++){
		gsl_matrix_set( work, k, l, gsl_ran_ugaussian(r) );
	}
}
gsl_matrix_memcpy(result,scale);
gsl_linalg_cholesky_decomp(result);
gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,result,work);
gsl_blas_dsyrk(CblasUpper,CblasNoTrans,1.0,work,0.0,result);

return 0;
}
