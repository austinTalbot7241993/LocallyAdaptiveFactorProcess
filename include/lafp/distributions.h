#ifndef MVNORM_H
#define MVNORM_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>

double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var);
int rmvt(const gsl_rng *r, const int n, const gsl_vector *location, const gsl_matrix *scale, const int dof, gsl_vector *result);
double dmvt(const int n, const gsl_vector *x, const gsl_vector *location, const gsl_matrix *scale, const int dof);
int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result);
int rwishart(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, gsl_matrix *result);

int ldmvnorm(int n,const gsl_vector *x,const gsl_vector *mean,const gsl_matrix *var,gsl_matrix *work,gsl_matrix *winv,double *ans);

#endif
