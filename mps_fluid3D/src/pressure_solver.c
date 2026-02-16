#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "pressure_solver.h"
#include "kernel.h"
#include "sim_config.h"

/* ベクトル内積 */
static double dot_product(const double *a, const double *b, int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
        sum += a[i] * b[i];
    return sum;
}

/* 行列-ベクトル積  y = A * x  (Aは n×n の1次元配列) */
static void mat_vec(const double *A, const double *x, double *y, int n)
{
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        y[i] = sum;
    }
}

/*
 * CG法 (共役勾配法)
 * A * x = b を解く (Aは対称正定値)
 * 戻り値: 反復回数
 */
static int solve_cg(const double *A, double *x, const double *b,
                    int n, int max_iter, double tol)
{
    double *r  = (double *)malloc(n * sizeof(double));
    double *p  = (double *)malloc(n * sizeof(double));
    double *Ap = (double *)malloc(n * sizeof(double));

    /* 初期値 x = 0, r = b, p = b */
    memset(x, 0, n * sizeof(double));
    memcpy(r, b, n * sizeof(double));
    memcpy(p, r, n * sizeof(double));

    double rr = dot_product(r, r, n);
    int iter;

    for (iter = 0; iter < max_iter; iter++) {
        if (sqrt(rr) < tol) break;

        mat_vec(A, p, Ap, n);

        double pAp = dot_product(p, Ap, n);
        if (fabs(pAp) < 1.0e-30) break;

        double alpha = rr / pAp;

        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rr_new = dot_product(r, r, n);
        double beta = rr_new / rr;

        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }

        rr = rr_new;
    }

    free(r);
    free(p);
    free(Ap);

    return iter;
}

/*
 * 不完全コレスキー分解 IC(0)
 * 対称行列 A に対して下三角行列 L を計算 (A ≈ L * L^T)
 * A の非零パターンのみに L の要素を生成する
 * L は n×n の密配列（下三角部分のみ使用）
 */
static void ic_factorize(const double *A, double *L, int n)
{
    memset(L, 0, (size_t)n * n * sizeof(double));

    for (int k = 0; k < n; k++) {
        /* 対角要素 L[k][k] = sqrt(A[k][k] - sum_{j<k} L[k][j]^2) */
        double sum = 0.0;
        for (int j = 0; j < k; j++) {
            if (fabs(A[k * n + j]) > 0.0) {
                sum += L[k * n + j] * L[k * n + j];
            }
        }
        double diag = A[k * n + k] - sum;
        if (diag <= 0.0) {
            /* 対角が非正の場合、元の対角値で代替 */
            diag = A[k * n + k];
        }
        L[k * n + k] = sqrt(diag);

        double Lkk_inv = 1.0 / L[k * n + k];

        /* 非対角要素 L[i][k] (i > k) */
        for (int i = k + 1; i < n; i++) {
            /* A[i][k] が 0 なら L[i][k] = 0 (IC(0)のスパースパターン保持) */
            if (fabs(A[i * n + k]) < 1.0e-30) continue;

            double s = 0.0;
            for (int j = 0; j < k; j++) {
                /* L[i][j] と L[k][j] が共に非零の場合のみ加算 */
                if (fabs(A[i * n + j]) > 0.0 && fabs(A[k * n + j]) > 0.0) {
                    s += L[i * n + j] * L[k * n + j];
                }
            }
            L[i * n + k] = (A[i * n + k] - s) * Lkk_inv;
        }
    }
}

/*
 * 前進代入  L * y = r
 * L は下三角行列 (n×n 密配列)
 */
static void forward_solve(const double *L, const double *r, double *y, int n)
{
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i * n + j] * y[j];
        }
        y[i] = (r[i] - sum) / L[i * n + i];
    }
}

/*
 * 後退代入  L^T * z = y
 * L は下三角行列 (n×n 密配列)
 */
static void backward_solve(const double *L, const double *y, double *z, int n)
{
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += L[j * n + i] * z[j];
        }
        z[i] = (y[i] - sum) / L[i * n + i];
    }
}

/*
 * 前処理の適用  M^{-1} * r = z  (M = L * L^T)
 * L * y = r → L^T * z = y
 */
static void precond_solve(const double *L, const double *r, double *z,
                          double *work, int n)
{
    forward_solve(L, r, work, n);
    backward_solve(L, work, z, n);
}

/*
 * ICCG法 (不完全コレスキー分解前処理付き共役勾配法)
 * A * x = b を解く (Aは対称正定値)
 * 戻り値: 反復回数
 */
static int solve_iccg(const double *A, double *x, const double *b,
                      int n, int max_iter, double tol)
{
    double *L    = (double *)calloc((size_t)n * n, sizeof(double));
    double *r    = (double *)malloc(n * sizeof(double));
    double *z    = (double *)malloc(n * sizeof(double));
    double *p    = (double *)malloc(n * sizeof(double));
    double *Ap   = (double *)malloc(n * sizeof(double));
    double *work = (double *)malloc(n * sizeof(double));

    if (!L || !r || !z || !p || !Ap || !work) {
        fprintf(stderr, "Error: ICCG solver memory allocation failed\n");
        free(L); free(r); free(z); free(p); free(Ap); free(work);
        memset(x, 0, n * sizeof(double));
        return 0;
    }

    /* 不完全コレスキー分解 */
    ic_factorize(A, L, n);

    /* 初期値 x = 0, r = b */
    memset(x, 0, n * sizeof(double));
    memcpy(r, b, n * sizeof(double));

    /* z = M^{-1} * r */
    precond_solve(L, r, z, work, n);

    /* p = z */
    memcpy(p, z, n * sizeof(double));

    double rz = dot_product(r, z, n);
    int iter;

    for (iter = 0; iter < max_iter; iter++) {
        double r_norm = sqrt(dot_product(r, r, n));
        if (r_norm < tol) break;

        mat_vec(A, p, Ap, n);

        double pAp = dot_product(p, Ap, n);
        if (fabs(pAp) < 1.0e-30) break;

        double alpha = rz / pAp;

        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        /* z = M^{-1} * r */
        precond_solve(L, r, z, work, n);

        double rz_new = dot_product(r, z, n);
        double beta = rz_new / rz;

        for (int i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }

        rz = rz_new;
    }

    free(L);
    free(r);
    free(z);
    free(p);
    free(Ap);
    free(work);

    return iter;
}

/*
 * 圧力ポアソン方程式を構築・求解
 */
void solve_pressure(ParticleSystem *ps, NeighborList *nl)
{
    int n = ps->num;
    double re = g_config->influence_radius;
    double n0 = ps->n0;
    double lambda = ps->lambda;
    double coeff = 2.0 * DIM / (n0 * lambda);
    double dt = g_config->dt;
    double dt2 = dt * dt;

    /* 未知数のマッピング: 内部流体粒子のみ */
    int *eq_idx = (int *)malloc(n * sizeof(int));
    int n_eq = 0;
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type == FLUID_PARTICLE && !ps->particles[i].on_surface) {
            eq_idx[i] = n_eq++;
        } else {
            eq_idx[i] = -1;
        }
    }

    if (n_eq == 0) {
        for (int i = 0; i < n; i++)
            ps->particles[i].pressure = 0.0;
        free(eq_idx);
        return;
    }

    /* 係数行列 M と右辺ベクトル c を構築 */
    double *M = (double *)calloc((size_t)n_eq * n_eq, sizeof(double));
    double *c = (double *)calloc(n_eq, sizeof(double));
    double *x = (double *)calloc(n_eq, sizeof(double));

    if (!M || !c || !x) {
        fprintf(stderr, "Error: pressure solver memory allocation failed (n_eq=%d)\n", n_eq);
        free(eq_idx);
        free(M);
        free(c);
        free(x);
        return;
    }

    for (int i = 0; i < n; i++) {
        if (eq_idx[i] < 0) continue;
        int ei = eq_idx[i];

        double sum_w = 0.0;
        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                double diff = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += diff * diff;
            }
            double r = sqrt(r2);
            double w = kernel_weight(r, re);

            if (eq_idx[j] >= 0) {
                int ej = eq_idx[j];
                M[ei * n_eq + ej] = -coeff * w;
            }
            sum_w += w;
        }

        /* 対角要素 */
        M[ei * n_eq + ei] = coeff * sum_w;

        /* 右辺 */
        c[ei] = (g_config->density / dt2) * (ps->particles[i].n - n0) / n0;

        /* 緩和係数を適用 */
        c[ei] *= g_config->relaxation_coeff;
    }

    /* ソルバーの選択 */
    int iters;
    if (g_config->solver_type == 1) {
        iters = solve_iccg(M, x, c, n_eq,
                           g_config->cg_max_iter, g_config->cg_tolerance);
    } else {
        iters = solve_cg(M, x, c, n_eq,
                         g_config->cg_max_iter, g_config->cg_tolerance);
    }
    (void)iters;

    /* 結果を粒子に反映 */
    for (int i = 0; i < n; i++) {
        if (eq_idx[i] >= 0) {
            ps->particles[i].pressure = x[eq_idx[i]];
        } else {
            ps->particles[i].pressure = 0.0;
        }
    }

    free(eq_idx);
    free(M);
    free(c);
    free(x);
}

void clamp_negative_pressure(ParticleSystem *ps)
{
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].pressure < 0.0)
            ps->particles[i].pressure = 0.0;
    }
}
