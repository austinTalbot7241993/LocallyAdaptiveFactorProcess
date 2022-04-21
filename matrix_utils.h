#ifndef _ABT_MATRIX_UTILS_H
#define _ABT_MATRIX_UTILS_H

#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_eigen.h>
#include "gmlib.h"

//extern int _gm_errNum;

//#define GMERR(ENUM) {_gm_errNum=ENUM;goto ERROR;}
//#define GMERRH(FUNCNAME,RET) \
//ERROR:\
//	fprintf(stderr,"Error %d on line %d in function (%s) in file '%s'\n",\
//	_gm_errNum,__LINE__,FUNCNAME,__FILE__);\
//	return( RET )

#define GSVA(NAME,N) gsl_vector *NAME = gsl_vector_alloc(N)
#define GSMA(N1,N2) gsl_matrix_alloc(N1,N2)

//#define GSMA(NAME,N1,N2) gsl_matrix *NAME = Austin_gsl_matrix_alloc(N1,N2)
//#define GSVA(NAME,N) gsl_vector *NAME = Austin_gsl_matrix_alloc(N);

int gsl_vector_slice(gsl_vector *v,int start,int end,gsl_vector *out);
int gsl_matrix_slice(gsl_matrix *m,int x1,int x2,int y1,int y2,gsl_matrix *out);

int reverseVector(gsl_vector *v,gsl_vector *rv);

int mpmTd2(gsl_matrix *m);

//Error checking
int matrixCheckSize(gsl_matrix *m,int x,int y);
int vectorCheckSize(gsl_vector *v,int x);
int cmatrixIsZero(const gsl_matrix *m);
int cvectorIsZero(const gsl_vector *v);
int matrixIsZero(gsl_matrix *m);
int vectorIsZero(gsl_vector *v);

//Redefining my own allocation methods. Combined
//with the macros above these are interchangeable

gsl_vector * D_gsl_vector_alloc(int n);
gsl_matrix * D_gsl_matrix_alloc(int n1,int n2);

//Nice printing of matrices and vectors
int printGSLVector(gsl_vector *v);
int printGSLVectorT(gsl_vector *v);
int printGSLMatrix(gsl_matrix *m);

int printGSLVectorR(gsl_vector *v);
int printGSLVectorRT(gsl_vector *v);
int printGSLMatrixR(gsl_matrix *m);

int printGSLVectorComplex(gsl_vector_complex *v);
int printGSLVectorComplexT(gsl_vector_complex *v);
int printGSLMatrixComplex(gsl_matrix_complex *m);

int printGSLVectorSize(gsl_vector *v);
int printGSLMatrixSize(gsl_matrix *m);

int meanVector(gsl_vector *v,double *mean);
int meanMatrix(gsl_matrix *m,int axis,gsl_vector *v);

int demeanVector(gsl_vector *v,gsl_vector *v1);
int demeanMatrix(gsl_matrix *m,int axis,gsl_matrix *m1);

//This just lets me allocate special types of matrices
gsl_matrix * gsl_matrix_set_diagonal(gsl_vector *v);
gsl_matrix * gsl_matrix_set_blockDiagonal2(gsl_matrix *m1,gsl_matrix *m2);
gsl_matrix * gsl_matrix_set_blockDiagonal3(gsl_matrix *m1,gsl_matrix *m2,gsl_matrix *m3);
gsl_vector * gsl_vector_set_ones(int N);
gsl_matrix * gsl_matrix_set_ones(int N1,int N2);

int getDiag(gsl_matrix *m,gsl_vector *v);
int getOffDiag(gsl_matrix *m,gsl_vector *v,int offset);

int setDiag(gsl_matrix *m,gsl_vector *v);
int setOffDiag(gsl_matrix *m,gsl_vector *v,int offset);
int setDiagConst(gsl_matrix *m,double d);
int setOffDiagConst(gsl_matrix *m,double d,int offset);


int sqrtVector(gsl_vector *v);
int sqrtMatrix(gsl_matrix *m);

int repmat(gsl_vector *v,int d1,int d2,gsl_matrix *m);
//Random number generator


int randnV(gsl_vector *v,gsl_rng *rand,double mean,double std);
int randnM(gsl_matrix *m,gsl_rng *rand,double mean,double std);

int randgV(gsl_vector *v,gsl_rng *rand,double rate,double scale);
int randgM(gsl_matrix *m,gsl_rng *rand,double rate,double scale);

int randigV(gsl_vector *v,gsl_rng *rand,double rate,double scale);
int randigM(gsl_matrix *m,gsl_rng *rand,double rate,double scale);

int randDiagN(gsl_matrix *m,gsl_rng *rand,double mean,double std);
int randDiagG(gsl_matrix *m,gsl_rng *rand,double rate,double scale);
int randDiagIG(gsl_matrix *m,gsl_rng *rand,double rate,double scale);

//Checking for characteristics

int positiveVector(gsl_vector *v);
int positiveMatrix(gsl_matrix *m);
int positiveDefinite(gsl_matrix *m);
int symmetricMatrix(gsl_matrix *m);

//Kronecker products
int kron_mm(gsl_matrix *m1,gsl_matrix *m2,gsl_matrix *result);
int kron_mv(gsl_matrix *m,gsl_vector *v,gsl_matrix *result);
int kron_vm(gsl_vector *v,gsl_matrix *m,gsl_matrix *result);
int kron_vv(gsl_vector *v1,gsl_vector *v2,gsl_matrix *reult);

int hankel(gsl_vector *c,gsl_vector *r,gsl_matrix *out);

//Sum columns and rows
int sum_col(gsl_matrix *m,gsl_vector *v);
int sum_row(gsl_matrix *m,gsl_vector *v);

int sum_absRow(gsl_matrix *m,gsl_vector *v);
int sum_absCol(gsl_matrix *m,gsl_vector *v);

int sum_abs2Row(gsl_matrix *m,gsl_vector *v);
int sum_abs2Col(gsl_matrix *m,gsl_vector *v);

int reshape(gsl_vector *v,int s1,int s2,gsl_matrix *m);

//These methods are miscellaneous functions
int cumProd(gsl_vector *v,gsl_vector *result);
int cumSum(gsl_vector *v,gsl_vector *result);

double vecProd(gsl_vector *v);
double vecSum(gsl_vector *v);


double trace(gsl_matrix *m);

//Difference of vectors
int diff(gsl_vector *v,gsl_vector *d);

int vector_diff2(gsl_vector *v1,gsl_vector *v2,double *result);

/*Here are the common dot products*/
int dot_vv(gsl_vector *v1,gsl_vector *v2,double* result);
int dot_mv(gsl_matrix *m,gsl_vector *v,gsl_vector *result);
int dot_vm(gsl_vector *v,gsl_matrix *m,gsl_vector *result);
int dot_mm(gsl_matrix *m1,gsl_matrix *m2,gsl_matrix *result);
int dot_amaT(gsl_matrix *a,gsl_matrix *m,gsl_matrix *result);
int dot_aTma(gsl_matrix *a,gsl_matrix *m,gsl_matrix *result);
int dot_vTmv(gsl_vector *v,gsl_matrix *m,double *result);
int dot_vvT(gsl_vector *v,gsl_matrix *result);
int dot_mTm(gsl_matrix *m,gsl_matrix *result);
int dot_mmT(gsl_matrix *m,gsl_matrix *result);

int vector_outer(gsl_vector *v,gsl_matrix *m);

/*Here we have a bunch of norms for vectors and matrices*/
int vnorm2(gsl_vector *v,double *result);
int vnorm1(gsl_vector *v,double *result);
int vnorm00(gsl_vector *v,double *result);
int vnormp(gsl_vector *v,double p,double *result);
int mnorm1(gsl_matrix *m,double *result);
int mnormF(gsl_matrix *m,double *result);
int mnorm00(gsl_matrix *m,double *result);
int conditionNumber(gsl_matrix *m,double *result);
/*Here are methods for solving systems of equations*/

/*Here are methods for getting matrix inverses
  transposes and hermetian stuff			*/

int inv_LU(gsl_matrix *m,gsl_matrix *minv);
int inv_chol(gsl_matrix *m,gsl_matrix *minv);
int solve_LU(gsl_matrix *m,gsl_vector *b,gsl_vector *result);
int svd_LU(gsl_matrix *m,gsl_matrix *U,gsl_matrix *V,gsl_vector *S);
int det_LU(gsl_matrix *A,double *result);
//Allocate memory

int gTranspose(gsl_matrix *m,gsl_matrix *mT);
int gHermetian(gsl_matrix_complex *c,gsl_matrix_complex *cH);


//eigenvalue methods. Starting with a t indicates
//that memory is allocated within. This means its 
//good for debugging. However, the ones without a
//t are better for high performance code
//Always sorted in descending value


int tgetSEigenvalues(gsl_matrix *m,gsl_vector *v);
int tgetSEigenvectors(gsl_matrix *m,gsl_vector *eigenvalues,gsl_matrix *eigenvectors);
int getSEigenvalues(gsl_matrix *m,gsl_matrix *mcpy,gsl_eigen_symm_workspace *w,gsl_vector *eval);
int getSEigenvectors(gsl_matrix *m,gsl_matrix *mcpy,gsl_eigen_symmv_workspace *w,gsl_vector *eval,gsl_matrix *evec);

int tgetEigenvalues(gsl_matrix *m,gsl_vector_complex *v);
int tgetEigenvectors(gsl_matrix *m,gsl_vector_complex *eigenvalues,gsl_matrix_complex *eigenvectors);
int getEigenvalues(gsl_matrix *m,gsl_matrix *mcpy,gsl_eigen_nonsymm_workspace *w,gsl_vector_complex *eval);
int getEigenvectors(gsl_matrix *m,gsl_matrix *mcpy,gsl_eigen_nonsymmv_workspace *w,gsl_vector_complex *eval,gsl_matrix_complex *evec);

int minEigenRadius(gsl_matrix *m,double *rad);
int maxEigenRadius(gsl_matrix *m,double *rad);

#define DOT_MV(m,v,r) 		gsl_blas_dgemv(CblasNoTrans,1.0,m,v,0.0,r);
#define DOT_VV(v1,v2,r)		gsl_blas_ddot(v1,v2,r);
#define DOT_MM(m1,m2,r)		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1,m2,0.0,r);

#endif
