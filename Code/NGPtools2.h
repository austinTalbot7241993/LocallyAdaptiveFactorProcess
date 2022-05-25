#ifndef _ABT_NGPTOOLS2_H
#define _ABT_NGPTOOLS2_H

#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include "mdarray.h"


int NGP_Z(int Np,int Ns,marray3d *Z);
int NGP_G(int Ns,double delta,int approx,gsl_matrix *G);
int NGP_W(int Ns,double delta,gsl_vector *varU,gsl_vector *varA,int approx,gsl_matrix *W);
int NGP_assemble_matrices(int Np,int Ns,gsl_vector *delta,double *sigEps,gsl_vector *sigU,gsl_vector *sigA,double *sigMu,double *sigAlph,int approx,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,gsl_vector *a0,gsl_matrix *P0);

#endif
