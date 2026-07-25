#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <lafp/gmlib.h>
#include <lafp/matrix_utils.h>
#include <lafp/mdarray_complex.h>
#include <lafp/KalmanFilter2c.h>

Kalman2c * Kalman2c_New(void) {
    Kalman2c *self = (Kalman2c *)GM_Malloc(sizeof(Kalman2c));
    if (self == NULL) GMERR(-11);
    if (Kalman2c_init(self)) GMERR(-21);
    return self;
GMERRH("Kalman2c_New", NULL);
}

int Kalman2c_init(Kalman2c *self) {
    if (self == NULL) GMERR(-11);
    memset(self, '\0', sizeof(Kalman2c));
    return 0;
GMERRH("Kalman2c_init", 1);
}

int Kalman2c_free(Kalman2c *s) {
    if (s == NULL) GMERR(-11);
    if (GM_FreezeListWithMethod(&GM_FreeGSLVector, "x",
                                &s->yt_Vp, &s->at_Vm, &s->zat_Vp, &s->TTat_Vm,
                                &s->KtVt_Vm, NULL)) GMERR(-21);
    if (GM_FreezeListWithMethod(&GM_FreeGSLMatrix, "x",
                                &s->a_Mtm, &s->RRQQ,
                                &s->ZPt_Mpm, &s->Zt_Mmp,
                                &s->ZPZ_Mpp, &s->Finvt_Mpp, &s->TTPt_Mmm,
                                &s->TTPtZt_Mmp, &s->Kt_Mmp, &s->KtZ_Mmm,
                                &s->L_Mmm, &s->Lt_Mmm, &s->TTPtLt_Mmm,
                                &s->RQR, &s->RRT, NULL)) GMERR(-31);
    if (marray3d_complex_free(s->P_Atmm)) GMERR(-41);

    GM_Free(s);
    return 0;
GMERRH("Kalman2c_free", 1);
}

#define V_ALLOC(V, SIZE) self->V = gsl_vector_complex_alloc(SIZE)
#define M_ALLOC(M, SIZE1, SIZE2) self->M = gsl_matrix_complex_alloc(SIZE1, SIZE2)
#define A_ALLOC3(A, SIZE1, SIZE2, SIZE3) self->A = marray3d_complex_alloc(SIZE1, SIZE2, SIZE3)

int Kalman2c_construct(Kalman2c *self, gsl_matrix_complex *y, gsl_vector_complex *a_init,
                       gsl_matrix_complex *P_init, marray3d_complex *Z, marray3d_complex *H,
                       marray3d_complex *T, marray3d_complex *R, marray3d_complex *Q,
                       gsl_matrix_complex *v, marray3d_complex *K, marray3d_complex *Finv) {
    if (self == NULL || y == NULL || a_init == NULL || Q == NULL) GMERR(-1);
    self->Nt = y->size1;
    self->Np = y->size2;
    self->Nm = a_init->size;
    self->Nr = Q->d2;

    int t = self->Nt;
    int p = self->Np;
    int m = self->Nm;
    int r = self->Nr;

    // Allocate vectors
    self->a_init = a_init;
    V_ALLOC(yt_Vp, p);
    V_ALLOC(at_Vm, m);
    V_ALLOC(zat_Vp, p);
    V_ALLOC(TTat_Vm, m);
    V_ALLOC(KtVt_Vm, m);

    // Allocate matrices
    self->P_init = P_init;
    self->v = v;
    self->y = y;

    M_ALLOC(a_Mtm, t, m);
    M_ALLOC(ZPt_Mpm, p, m);
    M_ALLOC(Zt_Mmp, m, p);
    M_ALLOC(ZPZ_Mpp, p, p);
    M_ALLOC(Finvt_Mpp, p, p);
    M_ALLOC(TTPtZt_Mmp, m, p);
    M_ALLOC(TTPt_Mmm, m, m);
    M_ALLOC(Kt_Mmp, m, p);
    M_ALLOC(L_Mmm, m, m);
    M_ALLOC(KtZ_Mmm, m, m);
    M_ALLOC(Lt_Mmm, m, m);
    M_ALLOC(TTPtLt_Mmm, m, m);
    M_ALLOC(RRQQ, m, r);
    M_ALLOC(RQR, m, m);
    M_ALLOC(RRT, m, r);

    // Allocate array3ds
    A_ALLOC3(P_Atmm, t, m, m);
    self->T = T;
    self->R = R;
    self->Q = Q;
    self->K = K;
    self->Finv = Finv;
    self->H = H;
    self->Z = Z;

    return 0;
GMERRH("Kalman2c_construct", 1);
}

#undef V_ALLOC
#undef M_ALLOC
#undef A_ALLOC3

static int _Kalman2c_checkSizes(Kalman2c *self) {
    if (self == NULL) GMERR(-1);
    int Np = self->Np, Nm = self->Nm, Nr = self->Nr, Nt = self->Nt;

    if (marray3d_complex_checkYZ(self->Z, Np, Nm)) GMERR(-11);
    if (marray3d_complex_checkYZ(self->H, Np, Np)) GMERR(-21);
    if (marray3d_complex_checkYZ(self->T, Nm, Nm)) GMERR(-31);
    if (marray3d_complex_checkYZ(self->R, Nm, Nr)) GMERR(-41);
    if (marray3d_complex_checkYZ(self->Q, Nr, Nr)) GMERR(-51);

    if (self->v == NULL || self->v->size1 != (size_t)Nt || self->v->size2 != (size_t)Np) GMERR(-61);
    if (self->P_init == NULL || self->P_init->size1 != (size_t)Nm || self->P_init->size2 != (size_t)Nm) GMERR(-71);

    if (marray3d_complex_checkX1orT(self->Z, Nt)) GMERR(-81);
    if (marray3d_complex_checkX1orT(self->H, Nt)) GMERR(-91);
    if (marray3d_complex_checkX1orT(self->T, Nt)) GMERR(-101);
    if (marray3d_complex_checkX1orT(self->R, Nt)) GMERR(-111);
    if (marray3d_complex_checkX1orT(self->Q, Nt)) GMERR(-121);

    return 0;
GMERRH("_Kalman2c_checkSizes", 1);
}

static int _Kalman2c_zeroOutputs(Kalman2c *self) {
    if (self == NULL) GMERR(-1);
    gsl_matrix_complex_set_zero(self->v);
    gsl_matrix_complex_set_zero(self->a_Mtm);
    if (marray3d_complex_set_zero(self->K)) GMERR(-11);
    if (marray3d_complex_set_zero(self->Finv)) GMERR(-21);
    if (marray3d_complex_set_zero(self->P_Atmm)) GMERR(-31);

    return 0;
GMERRH("_Kalman2c_zeroOutputs", 1);
}

int Kalman2c_Filter(gsl_matrix_complex *y, gsl_vector_complex *a_init, gsl_matrix_complex *P_init,
                    marray3d_complex *Z, marray3d_complex *H, marray3d_complex *T,
                    marray3d_complex *R, marray3d_complex *Q, gsl_matrix_complex *v,
                    marray3d_complex *K, marray3d_complex *Finv) {
    Kalman2c *self = NULL;
    if ((self = Kalman2c_New()) == NULL) GMERR(-11);
    if (Kalman2c_construct(self, y, a_init, P_init, Z, H, T, R, Q, v, K, Finv)) GMERR(-21);
    if (_Kalman2c_checkSizes(self)) GMERR(-31);
    if (Kalman2c_operations(self)) GMERR(-41);
    if (Kalman2c_free(self)) GMERR(-51);
    return 0;
GMERRH("Kalman2c_Filter", 1);
}

#define A_ROW(A, M, T_IDX) M = self->A->matrixList[T_IDX]
#define A_SET_ROW(A, M, T_IDX) marray3d_complex_set_X(self->A, self->M, T_IDX)
#define M_GET_ROW(M_MAT, V_VEC, T_IDX) gsl_matrix_complex_get_row(self->V_VEC, self->M_MAT, T_IDX)
#define M_SET_ROW(M_MAT, V_VEC, T_IDX) gsl_matrix_complex_set_row(self->M_MAT, T_IDX, self->V_VEC)
#define V_SUB(VV, UU) gsl_vector_complex_sub(self->VV, self->UU)
#define V_ADD(VV, UU) gsl_vector_complex_add(self->VV, self->UU)
#define M_ADD(MM, NN) gsl_matrix_complex_add(self->MM, self->NN)
#define M_SUB(MM, NN) gsl_matrix_complex_sub(self->MM, self->NN)

int Kalman2c_operations(Kalman2c *self) {
    int t;
    if (self == NULL) GMERR(-1);
    int Nt = self->Nt;

    int Ht = ((self->H->d1) > 1);
    int Tt = ((self->T->d1) > 1);
    int Rt = ((self->R->d1) > 1);
    int Qt = ((self->Q->d1) > 1);
    int Zt = ((self->Z->d1) > 1);

    gsl_matrix_complex *HH, *TT, *RR, *QQ, *ZZ, *Pt_Mmm;
    gsl_complex alpha_one = gsl_complex_rect(1.0, 0.0);
    gsl_complex beta_zero = gsl_complex_rect(0.0, 0.0);

    if (_Kalman2c_zeroOutputs(self)) GMERR(-11);

    if (M_SET_ROW(a_Mtm, a_init, 0)) GMERR(-21);
    if (A_SET_ROW(P_Atmm, P_init, 0)) GMERR(-31);
    A_ROW(H, HH, 0);
    A_ROW(T, TT, 0);
    A_ROW(R, RR, 0);
    A_ROW(Q, QQ, 0);
    A_ROW(Z, ZZ, 0);
    A_ROW(P_Atmm, Pt_Mmm, 0);
    gsl_matrix_complex_memcpy(self->Zt_Mmp, ZZ);

    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, RR, QQ, beta_zero, self->RRQQ);
    gsl_matrix_complex_memcpy(self->RRT, RR);

    for (t = 0; t < Nt; t++) {
        if (Ht) A_ROW(H, HH, t);
        if (Tt) A_ROW(T, TT, t);
        if (Rt) A_ROW(R, RR, t);
        if (Qt) A_ROW(Q, QQ, t);
        if (Zt) A_ROW(Z, ZZ, t);
        if (Zt) gsl_matrix_complex_transpose_memcpy(self->Zt_Mmp, ZZ);

        if (M_GET_ROW(a_Mtm, at_Vm, t)) GMERR(-141);
        if (M_GET_ROW(y, yt_Vp, t)) GMERR(-151);

        gsl_blas_zgemv(CblasNoTrans, alpha_one, ZZ, self->at_Vm, beta_zero, self->zat_Vp);
        if (V_SUB(yt_Vp, zat_Vp)) GMERR(-171);
        if (M_SET_ROW(v, yt_Vp, t)) GMERR(-181);

        A_ROW(P_Atmm, Pt_Mmm, t);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, ZZ, Pt_Mmm, beta_zero, self->ZPt_Mpm);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, self->ZPt_Mpm, self->Zt_Mmp, beta_zero, self->ZPZ_Mpp);
        if (gsl_matrix_complex_add(self->ZPZ_Mpp, HH)) GMERR(-221);

        if (A_SET_ROW(Finv, Finvt_Mpp, t)) GMERR(-241);

        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, TT, Pt_Mmm, beta_zero, self->TTPt_Mmm);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, self->TTPt_Mmm, self->Zt_Mmp, beta_zero, self->TTPtZt_Mmp);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, self->TTPtZt_Mmp, self->Finvt_Mpp, beta_zero, self->Kt_Mmp);
        if (A_SET_ROW(K, Kt_Mmp, t)) GMERR(-281);

        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, self->Kt_Mmp, ZZ, beta_zero, self->KtZ_Mmm);
        if (gsl_matrix_complex_memcpy(self->L_Mmm, TT)) GMERR(-301);
        if (M_SUB(L_Mmm, KtZ_Mmm)) GMERR(-311);

        if (t + 1 < Nt) {
            gsl_blas_zgemv(CblasNoTrans, alpha_one, TT, self->at_Vm, beta_zero, self->TTat_Vm);
            gsl_blas_zgemv(CblasNoTrans, alpha_one, self->Kt_Mmp, self->yt_Vp, beta_zero, self->KtVt_Vm);
            if (V_ADD(TTat_Vm, KtVt_Vm)) GMERR(-341);
            if (M_SET_ROW(a_Mtm, TTat_Vm, t + 1)) GMERR(-351);

            gsl_matrix_complex_transpose_memcpy(self->Lt_Mmm, self->L_Mmm);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, self->TTPt_Mmm, self->Lt_Mmm, beta_zero, self->TTPtLt_Mmm);
            if (Rt || Qt) gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, RR, QQ, beta_zero, self->RRQQ);
            if (Rt) gsl_matrix_complex_transpose_memcpy(self->RRT, RR);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, self->RRQQ, self->RRT, beta_zero, self->RQR);
            if (M_ADD(TTPtLt_Mmm, RQR)) GMERR(-411);
            if (A_SET_ROW(P_Atmm, TTPtLt_Mmm, t + 1)) GMERR(-421);
        }
    }
    return 0;
GMERRH("Kalman2c_operations", 1);
}

#undef A_ROW
#undef A_SET_ROW
#undef M_GET_ROW
#undef M_SET_ROW
#undef V_SUB
#undef V_ADD
#undef M_ADD
#undef M_SUB

int Kalman2c_writeV(Kalman2c *self, char *name) {
    FILE *fp;
    if (self == NULL) GMERR(-1);
    fp = fopen(name, "w");
    if (fp == NULL) GMERR(-3);
    if (gsl_matrix_complex_fprintf(fp, self->v, "%.15f")) GMERR(-11);
    fclose(fp);
    return 0;
GMERRH("Kalman2c_writeV", 1);
}

int Kalman2c_writeK(Kalman2c *self, char *name) {
    FILE *fp;
    if (self == NULL) GMERR(-1);
    fp = fopen(name, "w");
    if (fp == NULL) GMERR(-3);
    if (marray3d_complex_write(fp, self->K)) GMERR(-11);
    fclose(fp);
    return 0;
GMERRH("Kalman2c_writeK", 1);
}

int Kalman2c_writeFinv(Kalman2c *self, char *name) {
    FILE *fp;
    if (self == NULL) GMERR(-1);
    fp = fopen(name, "w");
    if (fp == NULL) GMERR(-3);
    if (marray3d_complex_write(fp, self->Finv)) GMERR(-11);
    fclose(fp);
    return 0;
GMERRH("Kalman2c_writeFinv", 1);
}

int Kalman2c_writeOutputs(Kalman2c *self, char *baseName) {
    char strV[16] = "_V.txt";
    char strK[16] = "_K.txt";
    char strFinv[16] = "_Finv.txt";
    char outV[1024], outK[1024], outFinv[1024];
    snprintf(outV, sizeof(outV), "%s%s", baseName, strV);
    snprintf(outK, sizeof(outK), "%s%s", baseName, strK);
    snprintf(outFinv, sizeof(outFinv), "%s%s", baseName, strFinv);

    if (self == NULL) GMERR(-1);
    if (Kalman2c_writeV(self, outV)) GMERR(-11);
    if (Kalman2c_writeK(self, outK)) GMERR(-21);
    if (Kalman2c_writeFinv(self, outFinv)) GMERR(-31);

    return 0;
GMERRH("Kalman2c_writeOutputs", 1);
}
