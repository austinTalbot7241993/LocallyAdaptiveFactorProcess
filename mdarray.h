#ifndef _ABT_MDARRAY_H
#define _ABT_MDARRAY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_permutation.h>

/*Code by Austin Bryan Talbot
 *        a.talbot@verizon.net
 *
 *Description: This code implements 3 and 4 dimensional
 *arrays. While there is such a thing as marray, I found
 *it quite unsatisfactory as it did not provide multiplication
 *or any of the other methods I cared about. In esscence, all
 *I needed was an organized list of matrices (my LAF code has
 *a covariance matrix for each time point). marray3d has 4 
 *atributes,d1,d2 and d3 which give the dimensions of the array,
 *and matrixList which contains the data. Note that the data is not
 *guaranteed to be contiguous since it is an array of pointers.
 *
 *Completed:12/24/2016
 *
 *Revision History:
 *
*/


typedef struct{
	size_t d1;
	size_t d2;
	size_t d3;
	gsl_matrix **matrixList;
} marray3d;

typedef struct{
	size_t d1;
	size_t d2;
	size_t d3;

} marray3dUpper;

typedef struct{
	size_t d1;
	size_t d2;
	size_t d3;
	size_t d4;
	marray3d **marrayList;

} marray4d;

typedef struct{
	size_t d1;
	size_t d2;
	size_t d3;
	size_t d4;

} marray4dUpper;

/***********************
* Methods for marray4d *
***********************/

//Allocation and free methods
marray4d * marray4d_alloc(size_t d1,size_t d2,size_t d3,size_t d4);
int marray4d_free(marray4d *self);

//Initialization methods
int marray4d_set_zero(marray4d *self);
int marray4d_set_constant(marray4d *self,double x);
int marray4d_set_identity(marray4d *self);
int marray4d_set_marray3d(marray4d *self,marray3d *m);

//Read and write methods
int marray4d_read(FILE *fp,marray4d *self);
int marray4d_write(FILE *fp,marray4d *self);
int marray4d_print(marray4d *self);
int marray4d_fancyWrite(FILE *fp,marray4d *self);

//Random initialization with a common random variable
int marray4d_set_randN(marray4d *self,double mu,double sigma2);
int marray4d_set_randG(marray4d *self,double rate,double scale);
int marray4d_set_randIG(marray4d *self,double rate,double scale);
int marray4d_set_randU(marray4d *self,double a,double b);

//Size comparisons
int marray4d_checkSizes(marray4d *self,size_t d1,size_t d2,size_t d3,size_t d4);
int marray4d_checkYZW(marray4d *self,size_t d2,size_t d3,size_t d4);
int marray4d_checkX1orT(marray4d *self,size_t t);
int marray4d_checkY1orT(marray4d *self,size_t t);
int marray4d_checkZ1orT(marray4d *self,size_t t);
int marray4d_checkW1orT(marray4d *self,size_t t);
int marray4d_compareSizes(marray4d *m1,marray4d *m2);

//Setters and getters
double marray4d_get(marray4d *self,size_t x,size_t y,size_t z,size_t w);
int marray4d_set(marray4d *self,size_t x,size_t y,size_t z,size_t w,double val);

int marray4d_set_ZW(marray4d *self,gsl_matrix *m,int x,int y);
int marray4d_set_X(marray4d *self,marray3d *m,int x);
int marray4d_get_ZW(marray4d *self,gsl_matrix *m,int x,int y);
int marray4d_get_X(marray4d *self,marray3d *m,int x);

/***********************
* Methods for marray3d *
***********************/

//Size comparisons
int marray3d_checkSizes(marray3d *self,size_t d1,size_t d2,size_t d3);
int marray3d_checkXY(marray3d *self,size_t d1,size_t d2);
int marray3d_checkYZ(marray3d *self,size_t d2,size_t d3);
int marray3d_checkXZ(marray3d *self,size_t d1,size_t d3);
int marray3d_checkX1orT(marray3d *self,size_t t);
int marray3d_checkY1orT(marray3d *self,size_t t);
int marray3d_checkZ1orT(marray3d *self,size_t t);
int marray3d_compareSizes(marray3d *m1,marray3d *m2);
int marray3d_compareXY(marray3d *m1,marray3d *m2);
int marray3d_compareYZ(marray3d *m1,marray3d *m2);
int marray3d_compareXZ(marray3d *m1,marray3d *m2);

//Allocation and free methods
marray3d * marray3d_alloc(size_t d1,size_t d2,size_t d3);
int marray3d_free(marray3d *self);

//Initialization methods
int marray3d_set_zero(marray3d *self);
int marray3d_set_constant(marray3d *self,double x);
int marray3d_set_identity(marray3d *self);
int marray3d_set_matrix(marray3d *self,gsl_matrix *m);

//Read and write methods
int marray3d_read(FILE *fp,marray3d *self);
int marray3d_write(FILE *fp,marray3d *self);
int marray3d_print(marray3d *self);
int marray3d_fancyWrite(FILE *fp,marray3d *self);

//Random initialization with a common random variable
int marray3d_set_randN(marray3d *self,gsl_rng *r,double mu,double sigma);
int marray3d_set_randG(marray3d *self,gsl_rng *r,double rate,double scale);
int marray3d_set_randIG(marray3d *self,gsl_rng *r,double rate,double scale);
int marray3d_set_randU(marray3d *self,gsl_rng *r,double a,double b);

//Random initializaiton with a function
//int marray3d_set_function(marray3d *self,fcnptr)

//Setters and getters
//Individual points
double marray3d_get(marray3d *self,size_t x,size_t y,size_t z);
int marray3d_set(marray3d *self,size_t x,size_t y,size_t z,double val);

//Slices
int marray3d_set_X(marray3d *self,gsl_matrix *m,int x);
int marray3d_set_Y(marray3d *self,gsl_matrix *m,int y);
int marray3d_set_Z(marray3d *self,gsl_matrix *m,int z);
int marray3d_get_X(marray3d *self,gsl_matrix *m,int x);
int marray3d_get_Y(marray3d *self,gsl_matrix *m,int y);
int marray3d_get_Z(marray3d *self,gsl_matrix *m,int z);

int marray3d_set_pencilX(marray3d *self,gsl_vector *v,int y,int z);
int marray3d_set_pencilY(marray3d *self,gsl_vector *v,int x,int z);
int marray3d_set_pencilZ(marray3d *self,gsl_vector *v,int x,int y);
int marray3d_get_pencilX(marray3d *self,gsl_vector *v,int y,int z);
int marray3d_get_pencilY(marray3d *self,gsl_vector *v,int x,int z);
int marray3d_get_pencilZ(marray3d *self,gsl_vector *v,int x,int y);

int marray3d_get_matrixSliceX(marray3d *self,gsl_matrix *m,int x,int y1,int y2,gsl_vector *work);
int marray3d_set_matrixSliceX(marray3d *self,gsl_matrix *m,int x,int y1,int y2,gsl_vector *work);

//Mathematics
int marray3d_LU(marray3d *self,marray3d *result,gsl_permutation *p,int *s);
int marray3d_cholesky(marray3d *self,marray3d *result);

//Math methods
int marray3d_add(marray3d *m1,marray3d *m2);
int marray3d_sub(marray3d *m1,marray3d *m2);
int marray3d_scale(marray3d *m1,double scalar);
int marray3d_cpy(marray3d *m1,marray3d *m2);
int marray3d_pointwiseMult(marray3d *m1,marray3d *m2);
int marray3d_pointwiseDivide(marray3d *m1,marray3d *m2);

int marray3d_mulX(marray3d *m1,marray3d *m2,marray3d *m3);
int marray3d_mulY(marray3d *m1,marray3d *m2,marray3d *m3);
int marray3d_mulZ(marray3d *m1,marray3d *m2,marray3d *m3);

int marray3d_mulMatX(marray3d *m1,gsl_matrix *m2,marray3d *m3);
int marray3d_mulMatY(marray3d *m1,gsl_matrix *m2,marray3d *m3);
int marray3d_mulMatZ(marray3d *m1,gsl_matrix *m2,marray3d *m3);


/*****************
* Hybrid methods *
*****************/

int marray34_compareYZW(marray4d *m4,marray3d *m3);

#endif
