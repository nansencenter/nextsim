/******************************************************************************
 *
 * File:        calcs.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "definitions.h"
#include "lapack.h"
#include "utils.h"
#include "calcs.h"

#define EPS 1.0e-8

static char doT = 'T';
static char noT = 'N';

/** Calculates trace of the product of two matrices.
 * A - n x m  (A[m][n])
 * B - m x n  (B[n][m])
 */
double traceprod(int transposeA, int transposeB, int m, int n, double** A, double** B)
{
    double trace = 0.0;
    int i, j;

#if 0
    if (!transposeA && !transposeB) {
        for (i = 0; i < n; ++i) {
            double* Bi = B[i];

            for (j = 0; j < m; ++j)
                trace += A[j][i] * Bi[j];
        }
    } else if (transposeA && transposeB) {
        for (i = 0; i < n; ++i) {
            double* Ai = A[i];

            for (j = 0; j < m; ++j)
                trace += Ai[j] * B[j][i];
        }
    } else if (transposeA || transposeB) {
        double* A0 = A[0];
        double* B0 = B[0];

        for (i = 0; i < n * m; ++i)
            trace += A0[i] * B0[i];
    }
#else
    if (!transposeA && !transposeB) {
        for (i = 0; i < n; ++i) {
            double* Aji = &A[0][i];
            double* Bij = &B[i][0];

            for (j = 0; j < m; ++j, Aji += n, Bij++)
                trace += *Aji * *Bij;
        }
    } else if (transposeA && transposeB) {
        for (i = 0; i < n; ++i) {
            double* Bji = &B[0][i];
            double* Aij = &A[i][0];

            for (j = 0; j < m; ++j, Bji += n, Aij++)
                trace += *Bji * *Aij;
        }
    } else if (transposeA || transposeB) {
        double* Aij = &A[0][0];
        double* Bij = &B[0][0];

        for (i = 0; i < n * m; ++i, Aij++, Bij++)
            trace += *Aij * *Bij;
    }
#endif

    return trace;
}

/** Calculates inverse of a symmetric matrix via Cholesky decomposition.
 * @param m - matrix size
 * @param S - input: matrix; output: S^-1
 * @return lapack_info from dgesvd_()
 */
static int invm(int m, double** S)
{
    char uplo = 'U';
    int lapack_info;
    int i, j;

    dpotrf_(&uplo, &m, S[0], &m, &lapack_info);
    if (lapack_info != 0)
        return lapack_info;
    dpotri_(&uplo, &m, S[0], &m, &lapack_info);
    if (lapack_info != 0)
        return lapack_info;
#if 0
    for (j = 1; j < m; ++j)
        for (i = 0; i < j; ++i)
            S[i][j] = S[j][i];
#else
    for (j = 1; j < m; ++j) {
        double* colj_rowi = S[j];
        double* rowj_coli = &S[0][j];

        for (i = 0; i < j; ++i, colj_rowi++, rowj_coli += m)
            *rowj_coli = *colj_rowi;
    }
#endif
    return 0;
}

/** Calculates inverse square root of a square matrix via SVD.
 * @param m - matrix size
 * @param S - input: matrix; output: S^-1
 * @return lapack_info from dgesvd_() (0 = success)
 */
int invsqrtm(int m, double** S)
{
    double** U = alloc2d(m, m, sizeof(double));
    double** Us = alloc2d(m, m, sizeof(double));
    double* sigmas = malloc(m * sizeof(double));
    int lwork = 10 * m;
    double* work = malloc(lwork * sizeof(double));
    char specU = 'A';           /* "all M columns of U are returned in array
                                 * * U" */
    char specV = 'N';           /* "no rows of V**T are computed" */
    int lapack_info;
    double a = 1.0;
    double b = 0.0;
    int i, j;

    dgesvd_(&specU, &specV, &m, &m, S[0], &m, sigmas, U[0], &m, NULL, &m, work, &lwork, &lapack_info);
    if (lapack_info != 0) {
        free(U);
        free(Us);
        free(sigmas);
        free(work);

        return lapack_info;
    }

    for (i = 0; i < m; ++i) {
        double* Ui = U[i];
        double* Usi = Us[i];
        double si_sqrt = sqrt(sigmas[i]);

        for (j = 0; j < m; ++j)
            Usi[j] = Ui[j] / si_sqrt;
    }
    dgemm_(&noT, &doT, &m, &m, &m, &a, Us[0], &m, U[0], &m, &b, S[0], &m);

    free(U);
    free(Us);
    free(sigmas);
    free(work);

    return 0;
}

/** Calculates M = inv(I + S' * S) [transpose = 0] or M = inv(I + S * S')
 ** [transpose = 1].
 * @param p - size(S, 1)
 * @param m - size(S, 2)
 * @param transpose - flag: whether to transpose S
 * @param S - input, p x m
 * @param M - output, m x m [transpose = 0], p x p [transpose = 1]
 * @return lapack_info (0 = success)
 */
static int calc_M(int p, int m, int transpose, double** S, double** M)
{
    int i;
    double a = 1.0;
    double b = 0.0;
    int lapack_info;

    assert(p > 0);

    if (!transpose) {
        /*
         * M = S' * S 
         */
        dgemm_(&doT, &noT, &m, &m, &p, &a, S[0], &p, S[0], &p, &b, M[0], &m);

        /*
         * M = I + S' * S
         */
        for (i = 0; i < m; ++i)
            M[i][i] += 1.0;

        /*
         * M = inv(I + S' * S) 
         */
        lapack_info = invm(m, M);
    } else {
        /*
         * M = S * S' 
         */
        dgemm_(&noT, &doT, &p, &p, &m, &a, S[0], &p, S[0], &p, &b, M[0], &p);
        /*
         * M = I + S * S'
         */
        for (i = 0; i < p; ++i)
            M[i][i] += 1.0;

        lapack_info = invm(p, M);
    }

    return lapack_info;
}

/** Calculates inverse _and_ inverse square root of a square matrix via SVD.
 * @param m - matrix size
 * @param S - input: matrix; output: S^-1
 * @param D - output: S^-1/2
 * @return lapack_info from dgesvd_() (0 = success)
 */
static int invsqrtm2(int m, double alpha, double** S, double** D)
{
    double** U = alloc2d(m, m, sizeof(double));
    double** Us1 = alloc2d(m, m, sizeof(double));
    double** Us2 = alloc2d(m, m, sizeof(double));
    double* sigmas = malloc(m * sizeof(double));
    int lwork = 10 * m;
    double* work = malloc(lwork * sizeof(double));
    char specU = 'A';           /* all M columns of U are returned in array U 
                                 */
    char specV = 'N';           /* no rows of V**T are computed */
    int lapack_info;
    double a = 1.0;
    double b = 0.0;
    int i, j;

    dgesvd_(&specU, &specV, &m, &m, S[0], &m, sigmas, U[0], &m, NULL, &m, work, &lwork, &lapack_info);
    if (lapack_info != 0) {
        free(U);
        free(Us1);
        free(Us2);
        free(sigmas);
        free(work);

        return lapack_info;
    }

    for (i = 0; i < m; ++i) {
        double* Ui = U[i];
        double* Us1i = Us1[i];
        double* Us2i = Us2[i];
        double si = sigmas[i];
        double si_sqrt = sqrt(1.0 - alpha + alpha * sigmas[i]);

        for (j = 0; j < m; ++j) {
            Us1i[j] = Ui[j] / si;
            Us2i[j] = Ui[j] / si_sqrt;
        }
    }
    dgemm_(&noT, &doT, &m, &m, &m, &a, Us1[0], &m, U[0], &m, &b, S[0], &m);
    dgemm_(&noT, &doT, &m, &m, &m, &a, Us2[0], &m, U[0], &m, &b, D[0], &m);

    free(U);
    free(Us1);
    free(Us2);
    free(sigmas);
    free(work);

    return 0;
}

/** Calculates G = inv(I + S' * S) * S' = S' * inv(I + S * S').
 */
void calc_G_denkf(int m, int p, double** S, int i, int j, double** G)
{
    double** M;
    double a = 1.0;
    double b = 0.0;
    int lapack_info;

    if (p < m) {
        M = alloc2d(p, p, sizeof(double));
        /*
         * M = inv(I + S * S')
         */
        lapack_info = calc_M(p, m, 1, S, M);
        if (lapack_info != 0)
            enkf_quit("dpotrf() or dpotri(): lapack_info = %d at (i, j) = (%d, %d)", lapack_info, i, j);

        /*
         * G = S' * inv(I + S * S') 
         */
        dgemm_(&doT, &noT, &m, &p, &p, &a, S[0], &p, M[0], &p, &b, G[0], &m);
    } else {
        /*
         * M = inv(I + S' * S)
         */
        M = alloc2d(m, m, sizeof(double));
        lapack_info = calc_M(p, m, 0, S, M);
        if (lapack_info != 0)
            enkf_quit("dpotrf() or dpotri(): lapack_info = %d at (i, j) = (%d, %d)", lapack_info, i, j);
        /*
         * G = inv(I + S * S') * S'
         */
        dgemm_(&noT, &doT, &m, &p, &m, &a, M[0], &m, S[0], &p, &b, G[0], &m);
    }

#if defined(CHECK_G)
    {
        int e, o;

        /*
         * check that columns of G sum up to 0
         */
        for (o = 1; o < p; ++o) {
            double* Go = G[o];
            double sumG = 0.0;
            double sumS = 0.0;

            for (e = 0; e < m; ++e) {
                sumG += Go[e];
                sumS += S[e][o];
            }

            if (fabs(sumG) > EPS)
                enkf_quit("inconsistency in G: column %d sums up to %.15f for (i, j) = (%d, %d); sum(S(%d,:) = %.15f)", o, sumG, i, j, o, sumS);
        }
    }
#endif

    free(M);
}

/** Calculates G = inv(I + S' * S) * S' and T = (I + S' * S)^-1/2.
 */
void calc_G_etkf(int m, int p, double** S, double alpha, int ii, int jj, double** G, double** T)
{
    double** M = alloc2d(m, m, sizeof(double));
    int lapack_info;
    int i;

    /*
     * dgemm stuff 
     */
    double a = 1.0;
    double b = 0.0;

    /*
     * M = S' * S 
     */
    dgemm_(&doT, &noT, &m, &m, &p, &a, S[0], &p, S[0], &p, &b, M[0], &m);

    /*
     * M = I + S' * S
     */
    for (i = 0; i < m; ++i)
        M[i][i] += 1.0;

    lapack_info = invsqrtm2(m, alpha, M, T);    /* M = M^-1, T = M^-1/2 */
    if (lapack_info != 0)
        enkf_quit("dgesvd(): lapack_info = %d at (i, j) = (%d, %d)", lapack_info, ii, jj);

    /*
     * G = inv(I + S * S') * S'
     */
    dgemm_(&noT, &doT, &m, &p, &m, &a, M[0], &m, S[0], &p, &b, G[0], &m);

#if defined(CHECK_G)
    /*
     * check that columns of G sum up to 0
     */
    {
        int e, o;

        for (o = 1; o < p; ++o) {
            double* Go = G[o];
            double sumG = 0.0;
            double sumS = 0.0;

            for (e = 0; e < m; ++e) {
                sumG += Go[e];
                sumS += S[e][o];
            }

            if (fabs(sumG) > EPS)
                enkf_quit("inconsistency in G: column %d sums up to %.15f for (i, j) = (%d, %d); sum(S(%d,:) = %.15f)", o, sumG, ii, jj, o, sumS);
        }
    }
#endif

    free(M);
}

/** Calculates X5 = G * s * 1' + T.
 * G is [m x p]
 * S is [p x m]
 * s is [p]
 * X5 is [m x m]
 */
void calc_X5_denkf(int m, int p, double** G, double** S, double* s, double alpha, int ii, int jj, double** X5)
{
    double* w = calloc(m, sizeof(double));
    double a, b;
    int i, j;

    /*
     * w = G * s 
     */
    if (enkf_nomeanupdate == 0) /* "normal" way */
        calc_w(m, p, G, s, w);

    /*
     * X5 = G * s * 1^T 
     */
    for (j = 0; j < m; ++j)
        memcpy(X5[j], w, m * sizeof(double));

    /*
     * X5 = G * s * 1^T - 1/2 * alpha * G * S 
     */
    a = -0.5 * alpha;
    b = 1.0;
    dgemm_(&noT, &noT, &m, &m, &p, &a, G[0], &m, S[0], &p, &b, X5[0], &m);

    /*
     * X5 = G * s * 1^T - 1/2 * alpha * G * S + I 
     */
    for (i = 0; i < m; ++i)
        X5[i][i] += 1.0;

#if defined(CHECK_X5)
    /*
     * check that columns of X5 sum up to 1 
     */
    for (i = 0; i < m; ++i) {
        double* X5i = X5[i];
        double sum = 0.0;

        for (j = 0; j < m; ++j)
            sum += X5i[j];

        if (fabs(sum - 1.0) > EPS)
            enkf_quit("inconsistency in X5: column %d sums up to %.15f for (i, j) = (%d, %d)", i, sum, ii, jj);
    }
#endif

    free(w);
}

/** Calculates X5 = G * s * 1' + T.
 * X5 = T on input
 */
void calc_X5_etkf(int m, int p, double** G, double* s, int ii, int jj, double** X5)
{
    double* w = calloc(m, sizeof(double));
    int i, j;

    /*
     * w = G * s 
     */
    if (enkf_nomeanupdate == 0) /* "normal" way */
        calc_w(m, p, G, s, w);

    /*
     * X5 = G * s * 1^T + T
     */
    for (i = 0; i < m; ++i) {
        double* X5i = X5[i];

        for (j = 0; j < m; ++j)
            X5i[j] += w[j];
    }

#if defined(CHECK_X5)
    /*
     * check that columns of X5 sum up to 1 
     */
    for (i = 0; i < m; ++i) {
        double* X5i = X5[i];
        double sum = 0.0;

        for (j = 0; j < m; ++j)
            sum += X5i[j];

        if (fabs(sum - 1.0) > EPS)
            enkf_quit("inconsistency in X5: column %d sums up to %.15f for (i, j) = (%d, %d)", i, sum, ii, jj);
    }
#endif

    free(w);
}

/** Calculates w = G * s.
 */
void calc_w(int m, int p, double** G, double* s, double* w)
{
    double a = 1.0;
    int inc = 1;
    double b = 0.0;

    dgemv_(&noT, &m, &p, &a, G[0], &m, s, &inc, &b, w, &inc);
}
