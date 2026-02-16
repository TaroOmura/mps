#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "pressure_solver.h"
#include "kernel.h"
#include "config.h"

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
 * 圧力ポアソン方程式を構築・求解
 *
 * ラプラシアンモデル:
 *   <∇²P>_i = (2d / (n0 * λ)) * Σ_{j≠i} (P_j - P_i) * w_ij
 *
 * ポアソン方程式:
 *   <∇²P>_i = -(ρ / Δt²) * (n*_i - n0) / n0
 *
 * 内部流体粒子のみを未知数として連立方程式を構築
 * 自由表面・壁粒子は P = 0 として扱う
 */
void solve_pressure(ParticleSystem *ps, NeighborList *nl)
{
    int n = ps->num;
    double re = INFLUENCE_RADIUS;
    double n0 = ps->n0;
    double lambda = ps->lambda;
    double coeff = 2.0 * DIM / (n0 * lambda);
    double dt2 = DT * DT;

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
        /* 全粒子が自由表面 or 壁 → 圧力はすべて0 */
        for (int i = 0; i < n; i++)
            ps->particles[i].pressure = 0.0;
        free(eq_idx);
        return;
    }

    /* 係数行列 M と右辺ベクトル c を構築 */
    /* 符号反転して M = -A (正定値行列) とする */
    double *M = (double *)calloc(n_eq * n_eq, sizeof(double));
    double *c = (double *)calloc(n_eq, sizeof(double));
    double *x = (double *)calloc(n_eq, sizeof(double));

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
                /* j は内部流体粒子 → 係数行列に登録 */
                int ej = eq_idx[j];
                M[ei * n_eq + ej] = -coeff * w;  /* -A の off-diagonal は負 */
            }
            /* j が境界粒子(P=0)の場合、RHSへの寄与は0 */
            sum_w += w;
        }

        /* 対角要素: -A_ii = coeff * Σ w_ij (正) */
        M[ei * n_eq + ei] = coeff * sum_w;

        /* 右辺: -b_i = (ρ / Δt²) * (n*_i - n0) / n0 */
        c[ei] = (DENSITY / dt2) * (ps->particles[i].n - n0) / n0;

        /* 緩和係数を適用 */
        c[ei] *= RELAXATION_COEFF;
    }

    /* CG法で求解 */
    int iters = solve_cg(M, x, c, n_eq, CG_MAX_ITER, CG_TOLERANCE);
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
