#ifndef _ABT_MDARRAY_COMPLEX_H
#define _ABT_MDARRAY_COMPLEX_H

/**
 * @file mdarray_complex.h
 * @brief 3D and 4D complex-valued multidimensional array abstractions using GSL.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_permutation.h>
#include <lafp/mdarray.h>

typedef struct {
    size_t d1;
    size_t d2;
    size_t d3;
    gsl_matrix_complex **matrixList;
} marray3d_complex;

typedef struct {
    size_t d1;
    size_t d2;
    size_t d3;
    size_t d4;
    marray3d_complex **marrayList;
} marray4d_complex;

/*******************************
* Methods for marray3d_complex *
*******************************/

int marray3d_complex_checkSizes(marray3d_complex *self, size_t d1, size_t d2, size_t d3);
int marray3d_complex_checkXY(marray3d_complex *self, size_t d1, size_t d2);
int marray3d_complex_checkYZ(marray3d_complex *self, size_t d2, size_t d3);
int marray3d_complex_checkXZ(marray3d_complex *self, size_t d1, size_t d3);
int marray3d_complex_checkX1orT(marray3d_complex *self, size_t t);
int marray3d_complex_checkY1orT(marray3d_complex *self, size_t t);
int marray3d_complex_checkZ1orT(marray3d_complex *self, size_t t);
int marray3d_complex_compareSizes(marray3d_complex *m1, marray3d_complex *m2);
int marray3d_complex_compareXY(marray3d_complex *m1, marray3d_complex *m2);
int marray3d_complex_compareYZ(marray3d_complex *m1, marray3d_complex *m2);
int marray3d_complex_compareXZ(marray3d_complex *m1, marray3d_complex *m2);

// Allocation and free methods
marray3d_complex * marray3d_complex_alloc(size_t d1, size_t d2, size_t d3);
int marray3d_complex_free(marray3d_complex *self);

// Initialization methods
int marray3d_complex_set_zero(marray3d_complex *self);
int marray3d_complex_set_constant(marray3d_complex *self, gsl_complex x);
int marray3d_complex_set_identity(marray3d_complex *self);
int marray3d_complex_set_matrix(marray3d_complex *self, gsl_matrix_complex *m);

// Read and write methods
int marray3d_complex_read(FILE *fp, marray3d_complex *self);
int marray3d_complex_write(FILE *fp, marray3d_complex *self);
int marray3d_complex_print(marray3d_complex *self);
int marray3d_complex_fancyWrite(FILE *fp, marray3d_complex *self);

// Setters and getters
gsl_complex marray3d_complex_get(marray3d_complex *self, size_t x, size_t y, size_t z);
int marray3d_complex_set(marray3d_complex *self, size_t x, size_t y, size_t z, gsl_complex val);

// Slices
int marray3d_complex_set_X(marray3d_complex *self, gsl_matrix_complex *m, int x);
int marray3d_complex_set_Y(marray3d_complex *self, gsl_matrix_complex *m, int y);
int marray3d_complex_set_Z(marray3d_complex *self, gsl_matrix_complex *m, int z);
int marray3d_complex_get_X(marray3d_complex *self, gsl_matrix_complex *m, int x);
int marray3d_complex_get_Y(marray3d_complex *self, gsl_matrix_complex *m, int y);
int marray3d_complex_get_Z(marray3d_complex *self, gsl_matrix_complex *m, int z);

int marray3d_complex_set_pencilX(marray3d_complex *self, gsl_vector_complex *v, int y, int z);
int marray3d_complex_set_pencilY(marray3d_complex *self, gsl_vector_complex *v, int x, int z);
int marray3d_complex_set_pencilZ(marray3d_complex *self, gsl_vector_complex *v, int x, int y);
int marray3d_complex_get_pencilX(marray3d_complex *self, gsl_vector_complex *v, int y, int z);
int marray3d_complex_get_pencilY(marray3d_complex *self, gsl_vector_complex *v, int x, int z);
int marray3d_complex_get_pencilZ(marray3d_complex *self, gsl_vector_complex *v, int x, int y);

int marray3d_complex_get_matrixSliceX(marray3d_complex *self, gsl_matrix_complex *m, int x, int y1, int y2, gsl_vector_complex *work);
int marray3d_complex_set_matrixSliceX(marray3d_complex *self, gsl_matrix_complex *m, int x, int y1, int y2, gsl_vector_complex *work);

// Mathematics
int marray3d_complex_LU(marray3d_complex *self, marray3d_complex *result, gsl_permutation *p, int *s);
int marray3d_complex_cholesky(marray3d_complex *self, marray3d_complex *result);

// Math methods
int marray3d_complex_add(marray3d_complex *m1, marray3d_complex *m2);
int marray3d_complex_sub(marray3d_complex *m1, marray3d_complex *m2);
int marray3d_complex_scale(marray3d_complex *m1, double scalar);
int marray3d_complex_cpy(marray3d_complex *m1, marray3d_complex *m2);
int marray3d_complex_pointwiseMult(marray3d_complex *m1, marray3d_complex *m2);
int marray3d_complex_pointwiseDivide(marray3d_complex *m1, marray3d_complex *m2);

int marray3d_complex_mulX(marray3d_complex *m1, marray3d_complex *m2, marray3d_complex *m3);
int marray3d_complex_mulY(marray3d_complex *m1, marray3d_complex *m2, marray3d_complex *m3);
int marray3d_complex_mulZ(marray3d_complex *m1, marray3d_complex *m2, marray3d_complex *m3);

int marray3d_complex_mulMatX(marray3d_complex *m1, gsl_matrix_complex *m2, marray3d_complex *m3);
int marray3d_complex_mulMatY(marray3d_complex *m1, gsl_matrix_complex *m2, marray3d_complex *m3);
int marray3d_complex_mulMatZ(marray3d_complex *m1, gsl_matrix_complex *m2, marray3d_complex *m3);

#endif
