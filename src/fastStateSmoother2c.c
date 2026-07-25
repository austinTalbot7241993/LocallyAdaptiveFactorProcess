#include <stdio.h>
#include <string.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>
#include <lafp/gmlib.h>
#include <lafp/matrix_utils.h>
#include <lafp/mdarray_complex.h>
#include <lafp/fastStateSmoother2c.h>

static int _FSS2c_checkSizes(FSS2c *self);

FSS2c * FSS2c_New(void) {
    FSS2c *self = (FSS2c *)GM_Malloc(sizeof(FSS2c));
    if (self == NULL) GMERR(-11);
    if (FSS2c_init(self)) GMERR(-21);
    return self;
GMERRH("FSS2c_New", NULL);
}

int FSS2c_init(FSS2c *self) {
    if (self == NULL) GMERR(-11);
    memset(self, '\0', sizeof(FSS2c));
    return 0;
GMERRH("FSS2c_init", 1);
}

int FSS2c_free(FSS2c *s) {
    if (s == NULL) GMERR(-11);
    if (GM_FreezeListWithMethod(&GM_FreeGSLVector, "x",
                                &s->vt_Vp, &s->Fvt_Vp, &s->Kttrt_Vp,
                                &s->rt_Vm, &s->r_init, &s->TTrt_Vm,
                                &s->ZZFvt_Vm, &s->Pr_Vm, &s->TTat_Vm,
                                &s->RQRrt_Vm, &s->at_Vm, NULL)) GMERR(-21);
    if (GM_FreezeListWithMethod(&GM_FreeGSLMatrix, "x",
                                &s->r, &s->ZZT_Mmp, &s->RRQQ, &s->RRt,
                                &s->RQR, &s->Ktt_Mpm,
                                &s->TTt_Mmm, NULL)) GMERR(-31);

    GM_Free(s);
    return 0;
GMERRH("FSS2c_free", 1);
}

#define M_ALLOC(M, SIZE1, SIZE2) self->M = gsl_matrix_complex_alloc(SIZE1, SIZE2)
#define V_ALLOC(V, SIZE) self->V = gsl_vector_complex_alloc(SIZE)

int FSS2c_construct(FSS2c *self, gsl_matrix_complex *v, marray3d_complex *K, marray3d_complex *Finv,
                    gsl_vector_complex *a_init, gsl_matrix_complex *P_init, marray3d_complex *Z,
                    marray3d_complex *T, marray3d_complex *R, marray3d_complex *Q, gsl_matrix_complex *alpha) {
    int t, p, m, r;
    if (self == NULL || v == NULL || a_init == NULL || Q == NULL) GMERR(-1);
    self->Nt = v->size1;
    self->Np = v->size2;
    self->Nm = a_init->size;
    self->Nr = Q->d2;

    t = self->Nt;
    p = self->Np;
    m = self->Nm;
    r = self->Nr;

    // Vectors
    self->a_init = a_init;

    V_ALLOC(vt_Vp, p);
    V_ALLOC(Fvt_Vp, p);
    V_ALLOC(Kttrt_Vp, p);
    V_ALLOC(rt_Vm, m);
    V_ALLOC(r_init, m);
    V_ALLOC(TTrt_Vm, m);
    V_ALLOC(ZZFvt_Vm, m);
    V_ALLOC(Pr_Vm, m);
    V_ALLOC(TTat_Vm, m);
    V_ALLOC(RQRrt_Vm, m);
    V_ALLOC(at_Vm, m);

    // Matrices
    self->v = v;
    self->P_init = P_init;
    self->alpha = alpha;

    M_ALLOC(r, t, m);
    M_ALLOC(ZZT_Mmp, m, p);
    M_ALLOC(RRQQ, m, r);
    M_ALLOC(RRt, r, m);
    M_ALLOC(RQR, m, m);
    M_ALLOC(Ktt_Mpm, p, m);
    M_ALLOC(TTt_Mmm, m, m);

    // Arrays
    self->Z = Z;
    self->K = K;
    self->Finv = Finv;
    self->T = T;
    self->R = R;
    self->Q = Q;

    return 0;
GMERRH("FSS2c_construct", 1);
}

#undef M_ALLOC
#undef V_ALLOC

static int _FSS2c_checkSizes(FSS2c *self) {
    int Np, Nm, Nr, Nt;
    if (self == NULL) GMERR(-1);
    Np = self->Np; Nm = self->Nm; Nr = self->Nr; Nt = self->Nt;

    if (marray3d_complex_checkSizes(self->K, Nt, Nm, Np)) GMERR(-11);
    if (marray3d_complex_checkSizes(self->Finv, Nt, Np, Np)) GMERR(-21);
    if (marray3d_complex_checkYZ(self->T, Nm, Nm)) GMERR(-31);
    if (marray3d_complex_checkYZ(self->Q, Nr, Nr)) GMERR(-41);
    if (marray3d_complex_checkYZ(self->Z, Np, Nm)) GMERR(-51);

    if (self->v == NULL || self->v->size1 != (size_t)Nt || self->v->size2 != (size_t)Np) GMERR(-61);
    if (self->P_init == NULL || self->P_init->size1 != (size_t)Nm || self->P_init->size2 != (size_t)Nm) GMERR(-71);
    if (self->alpha == NULL || self->alpha->size1 != (size_t)Nt || self->alpha->size2 != (size_t)Nm) GMERR(-81);
    if (self->a_init == NULL || self->a_init->size != (size_t)Nm) GMERR(-91);

    if (marray3d_complex_checkX1orT(self->T, Nt)) GMERR(-101);
    if (marray3d_complex_checkX1orT(self->R, Nt)) GMERR(-111);
    if (marray3d_complex_checkX1orT(self->Q, Nt)) GMERR(-121);
    if (marray3d_complex_checkX1orT(self->Z, Nt)) GMERR(-131);

    return 0;
GMERRH("_FSS2c_checkSizes", 1);
}

static int _FSS2c_zeroOutputs(FSS2c *self) {
    gsl_matrix_complex_set_zero(self->alpha);
    gsl_matrix_complex_set_zero(self->r);
    gsl_vector_complex_set_zero(self->r_init);
    return 0;
GMERRH("_FSS2c_zeroOutputs", 1);
}

int FSS2c_Smoother(gsl_matrix_complex *v, marray3d_complex *K, marray3d_complex *Finv,
                   gsl_vector_complex *a_init, gsl_matrix_complex *P_init,
                   marray3d_complex *Z, marray3d_complex *T, marray3d_complex *R,
                   marray3d_complex *Q, gsl_matrix_complex *alpha) {
    FSS2c *self = NULL;
    if ((self = FSS2c_New()) == NULL) GMERR(-11);
    if (FSS2c_construct(self, v, K, Finv, a_init, P_init, Z, T, R, Q, alpha)) GMERR(-21);
    if (_FSS2c_checkSizes(self)) GMERR(-31);
    if (FSS2c_operations(self)) GMERR(-41);
    if (FSS2c_free(self)) GMERR(-51);
    return 0;
GMERRH("FSS2c_Smoother", 1);
}

#define A_ROW(A, M, T_IDX) M = self->A->matrixList[T_IDX]
#define M_GET_ROW(M_MAT, V_VEC, T_IDX) gsl_matrix_complex_get_row(self->V_VEC, self->M_MAT, T_IDX)
#define M_SET_ROW(M_MAT, V_VEC, T_IDX) gsl_matrix_complex_set_row(self->M_MAT, T_IDX, self->V_VEC)
#define V_SUB(VV, UU) gsl_vector_complex_sub(self->VV, self->UU)
#define V_ADD(VV, UU) gsl_vector_complex_add(self->VV, self->UU)

int FSS2c_operations(FSS2c *self) {
    int t;
    if (self == NULL) GMERR(-1);
    int Nt = self->Nt;
    int Tt = ((self->T->d1) > 1);
    int Rt = ((self->R->d1) > 1);
    int Qt = ((self->Q->d1) > 1);
    int Zt = ((self->Z->d1) > 1);

    gsl_matrix_complex *ZZ, *TT, *RR, *QQ, *Kt_Mmp, *Finvt_Mpp;
    gsl_complex alpha_one = gsl_complex_rect(1.0, 0.0);
    gsl_complex beta_zero = gsl_complex_rect(0.0, 0.0);

    if (_FSS2c_zeroOutputs(self)) GMERR(-1);

    A_ROW(Z, ZZ, 0);
    A_ROW(T, TT, 0);
    A_ROW(R, RR, 0);
    A_ROW(Q, QQ, 0);

    gsl_matrix_complex_transpose_memcpy(self->ZZT_Mmp, ZZ);

    for (t = Nt - 1; t > -1; t--) {
        if (Tt) A_ROW(T, TT, t);
        if (Rt) A_ROW(R, RR, t);
        if (Qt) A_ROW(Q, QQ, t);
        if (Zt) A_ROW(Z, ZZ, t);
        if (Zt) gsl_matrix_complex_transpose_memcpy(self->ZZT_Mmp, ZZ);

        A_ROW(Finv, Finvt_Mpp, t);
        A_ROW(K, Kt_Mmp, t);
        if (M_GET_ROW(r, rt_Vm, t)) GMERR(-101);
        if (M_GET_ROW(v, vt_Vp, t)) GMERR(-111);

        if (gsl_blas_zgemv(CblasNoTrans, alpha_one, Finvt_Mpp, self->vt_Vp, beta_zero, self->Fvt_Vp)) GMERR(-121);
        gsl_matrix_complex_transpose_memcpy(self->Ktt_Mpm, Kt_Mmp);
        if (gsl_blas_zgemv(CblasNoTrans, alpha_one, self->Ktt_Mpm, self->rt_Vm, beta_zero, self->Kttrt_Vp)) GMERR(-141);
        if (V_SUB(Fvt_Vp, Kttrt_Vp)) GMERR(-151);

        if (gsl_blas_zgemv(CblasNoTrans, alpha_one, self->ZZT_Mmp, self->Fvt_Vp, beta_zero, self->ZZFvt_Vm)) GMERR(-161);
        gsl_matrix_complex_transpose_memcpy(self->TTt_Mmm, TT);

        if (gsl_blas_zgemv(CblasNoTrans, alpha_one, self->TTt_Mmm, self->rt_Vm, beta_zero, self->TTrt_Vm)) GMERR(-181);
        if (V_ADD(ZZFvt_Vm, TTrt_Vm)) GMERR(-191);

        if (t > 0) {
            if (M_SET_ROW(r, ZZFvt_Vm, t - 1)) GMERR(-201);
        } else {
            if (gsl_vector_complex_memcpy(self->r_init, self->ZZFvt_Vm)) GMERR(-211);
        }
    }

    if (gsl_blas_zgemv(CblasNoTrans, alpha_one, self->P_init, self->r_init, beta_zero, self->Pr_Vm)) GMERR(-221);
    if (V_ADD(Pr_Vm, a_init)) GMERR(-231);
    if (M_SET_ROW(alpha, Pr_Vm, 0)) GMERR(-241);

    if (gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, RR, QQ, beta_zero, self->RRQQ)) GMERR(-251);
    gsl_matrix_complex_transpose_memcpy(self->RRt, RR);
    if (gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, self->RRQQ, self->RRt, beta_zero, self->RQR)) GMERR(-271);

    for (t = 0; t < Nt - 1; t++) {
        if (Tt) A_ROW(T, TT, t);
        if (Rt) A_ROW(R, RR, t);
        if (Qt) A_ROW(Q, QQ, t);

        if (Rt || Qt) {
            if (gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, RR, QQ, beta_zero, self->RRQQ)) GMERR(-451);
            gsl_matrix_complex_transpose_memcpy(self->RRt, RR);
            if (gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha_one, self->RRQQ, self->RRt, beta_zero, self->RQR)) GMERR(-471);
        }

        if (M_GET_ROW(alpha, at_Vm, t)) GMERR(-281);
        if (M_GET_ROW(r, rt_Vm, t)) GMERR(-291);
        if (gsl_blas_zgemv(CblasNoTrans, alpha_one, TT, self->at_Vm, beta_zero, self->TTat_Vm)) GMERR(-301);
        if (gsl_blas_zgemv(CblasNoTrans, alpha_one, self->RQR, self->rt_Vm, beta_zero, self->RQRrt_Vm)) GMERR(-311);
        if (V_ADD(TTat_Vm, RQRrt_Vm)) GMERR(-321);
        if (M_SET_ROW(alpha, TTat_Vm, (t + 1))) GMERR(-331);
    }
    return 0;
GMERRH("FSS2c_operations", 1);
}

#undef A_ROW
#undef M_GET_ROW
#undef M_SET_ROW
#undef V_ADD
#undef V_SUB

int FSS2c_writeAlpha(FSS2c *self, char *name) {
    FILE *fp;
    if (self == NULL) GMERR(-1);
    fp = fopen(name, "w");
    if (fp == NULL) GMERR(-3);
    if (gsl_matrix_complex_fprintf(fp, self->alpha, "%.15f")) GMERR(-11);
    fclose(fp);
    return 0;
GMERRH("FSS2c_writeAlpha", 1);
}

int FSS2c_writeOutputs(FSS2c *self, char *baseName) {
    char strA[16] = "_Alpha.txt";
    char outA[1024];
    snprintf(outA, sizeof(outA), "%s%s", baseName, strA);
    if (self == NULL) GMERR(-1);
    if (FSS2c_writeAlpha(self, outA)) GMERR(-11);
    return 0;
GMERRH("FSS2c_writeOutputs", 1);
}
