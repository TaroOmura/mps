/* OpenMP parallelization added */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "pressure_solver.h"
#include "kernel.h"
#include "sim_config.h"

/* ------------------------------------------------------------------ */
/* CSR疎行列                                                            */
/* ------------------------------------------------------------------ */

typedef struct {
    int     n;        /* 行数 (= 列数) */
    int     nnz;      /* 非ゼロ要素数 */
    int    *row_ptr;  /* [n+1] */
    int    *col_idx;  /* [nnz], 各行内でソート済み */
    double *val;      /* [nnz] */
    double *diag;     /* [n] 対角要素 (IC用) */
} CsrMatrix;

static CsrMatrix *csr_alloc(int n, int nnz)
{
    CsrMatrix *A = (CsrMatrix *)malloc(sizeof(CsrMatrix));
    if (!A) return NULL;
    A->n       = n;
    A->nnz     = nnz;
    A->row_ptr = (int    *)malloc((n + 1) * sizeof(int));
    A->col_idx = (int    *)malloc(nnz     * sizeof(int));
    A->val     = (double *)malloc(nnz     * sizeof(double));
    A->diag    = (double *)calloc(n,         sizeof(double));
    if (!A->row_ptr || !A->col_idx || !A->val || !A->diag) {
        free(A->row_ptr); free(A->col_idx); free(A->val); free(A->diag);
        free(A);
        return NULL;
    }
    return A;
}

static void csr_free(CsrMatrix *A)
{
    if (!A) return;
    free(A->row_ptr); free(A->col_idx); free(A->val); free(A->diag);
    free(A);
}

/*
 * y = A * x (行列ベクトル積)
 *
 * OpenMP並列化: 外ループ（行i）を並列化する。
 * 各行の内ループ（列方向sum）は各スレッドが独立して逐次集計するため、
 * 加算順序は逐次版と完全に同じ → y[i] の値はビット同一。
 * データ競合なし: 各スレッドは独立した y[i] に書き込む。
 */
static void csr_mat_vec(const CsrMatrix *A, const double *x, double *y)
{
    int n = A->n;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int p = A->row_ptr[i]; p < A->row_ptr[i + 1]; p++)
            s += A->val[p] * x[A->col_idx[p]];
        y[i] = s;
    }
}

/*
 * ドット積 (内積) — 逐次実装
 *
 * 並列 reduction(+:s) は浮動小数点の加算順序を変えるため、
 * ICCGの alpha/beta 計算値がわずかに変わり、数百ステップ後に
 * 累積誤差が圧力発散を引き起こすことが確認された。
 * 逐次実装により逐次コードとビット同一の収束挙動を保証する。
 * (n_eq ≈ 13000 要素の逐次sumは ~3μs であり性能影響は無視できる)
 */
static double dot(const double *a, const double *b, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

/* ------------------------------------------------------------------ */
/* CG法 (疎行列版, 並列化)                                              */
/* ------------------------------------------------------------------ */

static int solve_cg_sparse(const CsrMatrix *A, double *x, const double *b,
                           int max_iter, double tol)
{
    int n = A->n;
    double *r  = (double *)malloc(n * sizeof(double));
    double *p  = (double *)malloc(n * sizeof(double));
    double *Ap = (double *)malloc(n * sizeof(double));
    if (!r || !p || !Ap) { free(r); free(p); free(Ap); return 0; }

    memset(x, 0, n * sizeof(double));
    memcpy(r, b, n * sizeof(double));
    memcpy(p, r, n * sizeof(double));

    double rr = dot(r, r, n);
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        if (sqrt(rr) < tol) break;
        csr_mat_vec(A, p, Ap);
        double pAp = dot(p, Ap, n);
        if (fabs(pAp) < 1.0e-30) break;
        double alpha = rr / pAp;

        /* x += alpha * p;  r -= alpha * Ap; (並列) */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rr_new = dot(r, r, n);
        if (rr < 1.0e-60) break;
        double beta = rr_new / rr;

        /* p = r + beta * p (並列) */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) p[i] = r[i] + beta * p[i];

        rr = rr_new;
    }
    free(r); free(p); free(Ap);
    return iter;
}

/* ------------------------------------------------------------------ */
/* IC(0)前処理 (疎行列CSR版)                                            */
/*                                                                      */
/* IC分解は本質的に逐次的な依存関係を持つ（前の行の結果が後の行に影響）。 */
/* そのため IC 分解の内部は逐次のまま維持する。                          */
/* 前進/後退代入も依存関係があり逐次処理。                               */
/* ------------------------------------------------------------------ */

typedef struct {
    int     n;
    int    *row_ptr;
    int    *col_idx;
    double *val;
    double *diag_inv;
} IcFactor;

static void ic_free(IcFactor *L)
{
    if (!L) return;
    free(L->row_ptr); free(L->col_idx); free(L->val); free(L->diag_inv);
    free(L);
}

static IcFactor *ic_factorize_sparse(const CsrMatrix *A)
{
    int n = A->n;

    int lt_nnz = 0;
    for (int i = 0; i < n; i++)
        for (int p = A->row_ptr[i]; p < A->row_ptr[i + 1]; p++)
            if (A->col_idx[p] < i) lt_nnz++;

    IcFactor *L = (IcFactor *)malloc(sizeof(IcFactor));
    if (!L) return NULL;
    L->n        = n;
    L->row_ptr  = (int    *)malloc((n + 1) * sizeof(int));
    L->col_idx  = (int    *)malloc((lt_nnz + 1) * sizeof(int));
    L->val      = (double *)calloc((lt_nnz + 1),  sizeof(double));
    L->diag_inv = (double *)malloc(n              * sizeof(double));
    if (!L->row_ptr || !L->col_idx || !L->val || !L->diag_inv) {
        ic_free(L); return NULL;
    }

    L->row_ptr[0] = 0;
    for (int i = 0; i < n; i++) {
        int cnt = 0;
        for (int p = A->row_ptr[i]; p < A->row_ptr[i + 1]; p++)
            if (A->col_idx[p] < i) cnt++;
        L->row_ptr[i + 1] = L->row_ptr[i] + cnt;
    }
    for (int i = 0; i < n; i++) {
        int pos = L->row_ptr[i];
        for (int p = A->row_ptr[i]; p < A->row_ptr[i + 1]; p++) {
            int j = A->col_idx[p];
            if (j < i) {
                L->col_idx[pos] = j;
                L->val[pos]     = A->val[p];
                pos++;
            }
        }
    }

    /* IC(0) 分解 (逐次: 前の行への依存があるため並列化不可) */
    for (int k = 0; k < n; k++) {
        double d = A->diag[k];
        for (int p = L->row_ptr[k]; p < L->row_ptr[k + 1]; p++)
            d -= L->val[p] * L->val[p];
        if (d <= 0.0) d = (A->diag[k] > 0.0) ? A->diag[k] : 1.0e-30;
        L->diag_inv[k] = 1.0 / sqrt(d);

        for (int pk = A->row_ptr[k]; pk < A->row_ptr[k + 1]; pk++) {
            int i = A->col_idx[pk];
            if (i <= k) continue;

            double s = A->val[pk];

            int pi = L->row_ptr[i];
            int pk2 = L->row_ptr[k];
            int pi_end  = L->row_ptr[i + 1];
            int pk2_end = L->row_ptr[k + 1];
            while (pi < pi_end && pk2 < pk2_end) {
                int ji = L->col_idx[pi];
                int jk = L->col_idx[pk2];
                if      (ji < jk) { pi++; }
                else if (ji > jk) { pk2++; }
                else { s -= L->val[pi] * L->val[pk2]; pi++; pk2++; }
            }

            for (int p = L->row_ptr[i]; p < L->row_ptr[i + 1]; p++) {
                if (L->col_idx[p] == k) {
                    L->val[p] = s * L->diag_inv[k];
                    break;
                }
                if (L->col_idx[p] > k) break;
            }
        }
    }

    return L;
}

/* 前進代入 L * y = r (逐次: データ依存あり) */
static void forward_solve_sparse(const IcFactor *L, const double *r, double *y)
{
    for (int i = 0; i < L->n; i++) {
        double s = r[i];
        for (int p = L->row_ptr[i]; p < L->row_ptr[i + 1]; p++)
            s -= L->val[p] * y[L->col_idx[p]];
        y[i] = s * L->diag_inv[i];
    }
}

/* 後退代入 L^T * z = y (逐次: データ依存あり) */
static void backward_solve_sparse(const IcFactor *L, const double *y, double *z)
{
    int n = L->n;
    for (int i = 0; i < n; i++) z[i] = y[i];
    for (int i = n - 1; i >= 0; i--) {
        z[i] *= L->diag_inv[i];
        for (int p = L->row_ptr[i]; p < L->row_ptr[i + 1]; p++)
            z[L->col_idx[p]] -= L->val[p] * z[i];
    }
}

/* ------------------------------------------------------------------ */
/* ICCG法 (疎行列版, 部分並列化)                                        */
/* ------------------------------------------------------------------ */

static int solve_iccg_sparse(const CsrMatrix *A, double *x, const double *b,
                             int max_iter, double tol)
{
    int n = A->n;
    IcFactor *L = ic_factorize_sparse(A);
    if (!L) return solve_cg_sparse(A, x, b, max_iter, tol);

    double *r  = (double *)malloc(n * sizeof(double));
    double *z  = (double *)malloc(n * sizeof(double));
    double *p  = (double *)malloc(n * sizeof(double));
    double *Ap = (double *)malloc(n * sizeof(double));
    double *tmp= (double *)malloc(n * sizeof(double));
    if (!r || !z || !p || !Ap || !tmp) {
        free(r); free(z); free(p); free(Ap); free(tmp);
        ic_free(L);
        return solve_cg_sparse(A, x, b, max_iter, tol);
    }

    memset(x, 0, n * sizeof(double));
    memcpy(r, b, n * sizeof(double));

    /* z = M^{-1} r,  p = z (逐次: 前進/後退代入) */
    forward_solve_sparse(L, r, tmp);
    backward_solve_sparse(L, tmp, z);
    memcpy(p, z, n * sizeof(double));

    double rz = dot(r, z, n);
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        if (sqrt(dot(r, r, n)) < tol) break;

        csr_mat_vec(A, p, Ap);  /* 並列 */
        double pAp = dot(p, Ap, n);  /* 並列 (reduction) */
        if (fabs(pAp) < 1.0e-30) break;
        double alpha = rz / pAp;

        /* x += alpha * p;  r -= alpha * Ap; (並列) */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        /* z = M^{-1} r (逐次) */
        forward_solve_sparse(L, r, tmp);
        backward_solve_sparse(L, tmp, z);

        double rz_new = dot(r, z, n);  /* 並列 (reduction) */
        if (fabs(rz) < 1.0e-60) break;
        double beta = rz_new / rz;

        /* p = z + beta * p (並列) */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) p[i] = z[i] + beta * p[i];

        rz = rz_new;
    }

    free(r); free(z); free(p); free(Ap); free(tmp);
    ic_free(L);
    return iter;
}

/* ------------------------------------------------------------------ */
/* ソース項                                                              */
/* ------------------------------------------------------------------ */

static double source_term_standard(int i, const ParticleSystem *ps,
                                   double n0, double dt, const double *vel_div)
{
    double dt2 = dt * dt;
    if (g_config->ppe_type == 1) {
        return -g_config->density / dt * vel_div[i]
               + g_config->gamma_ppe * (g_config->density / dt2)
                 * (ps->particles[i].n - n0) / n0;
    } else {
        return g_config->relaxation_coeff
               * (g_config->density / dt2) * (ps->particles[i].n - n0) / n0;
    }
}

static double source_term_hs(int i, const ParticleSystem *ps,
                             NeighborList *nl, double n0, double dt, double re)
{
    int max_nb = nl->max_neighbors;
    double sum_hst = 0.0;
    for (int k = 0; k < nl->count[i]; k++) {
        int j = nl->neighbors[i * max_nb + k];
        if (ps->particles[j].type == DUMMY_PARTICLE) continue;
        double dr[DIM], r2 = 0.0;
        for (int d = 0; d < DIM; d++) {
            dr[d] = ps->particles[j].pos[d] - ps->particles[i].pos[d];
            r2 += dr[d] * dr[d];
        }
        if (r2 < 1.0e-24) continue;
        double r = sqrt(r2), dv_dr = 0.0;
        for (int d = 0; d < DIM; d++)
            dv_dr += dr[d] * (ps->particles[j].vel[d] - ps->particles[i].vel[d]);
        sum_hst += re / (r * r * r) * dv_dr;
    }
    return -(g_config->density / n0 / dt) * sum_hst;
}

/* ------------------------------------------------------------------ */
/* 圧力ポアソン方程式の構築・求解                                          */
/* ------------------------------------------------------------------ */

void solve_pressure(ParticleSystem *ps, NeighborList *nl)
{
    int    n      = ps->num;
    double re     = g_config->influence_radius_lap;
    double n0     = ps->n0;
    double lambda = ps->lambda;
    double coeff  = 2.0 * DIM / (n0 * lambda);
    double dt     = g_config->dt;
    int    max_nb = nl->max_neighbors;

    /* Natsui型: 速度発散の事前計算
     * 各粒子iは独立してvel_div[i]を計算 → 並列化可能。
     */
    double *vel_div = NULL;
    if (g_config->ppe_type == 1 && g_config->hs_mode == 0) {
        vel_div = (double *)calloc(n, sizeof(double));
        if (!vel_div) {
            fprintf(stderr, "Error: velocity divergence buffer allocation failed\n");
            return;
        }

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) {
            if (ps->particles[i].type != FLUID_PARTICLE) continue;
            double div_u = 0.0;
            for (int k = 0; k < nl->count[i]; k++) {
                int j = nl->neighbors[i * max_nb + k];
                if (ps->particles[j].type == DUMMY_PARTICLE) continue;
                double r2 = 0.0, dr[DIM];
                for (int d = 0; d < DIM; d++) {
                    dr[d] = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                    r2 += dr[d] * dr[d];
                }
                if (r2 < 1.0e-24) continue;
                double dv_dr = 0.0;
                for (int d = 0; d < DIM; d++)
                    dv_dr += (ps->particles[j].vel[d] - ps->particles[i].vel[d]) * dr[d];
                div_u += dv_dr / r2 * kernel_weight(sqrt(r2), re);
            }
            vel_div[i] = (double)DIM / n0 * div_u;
        }
    }

    /* 未知数マッピング (逐次: prefix-sum的な n_eq カウントが必要) */
    int *eq_idx = (int *)malloc(n * sizeof(int));
    if (!eq_idx) { free(vel_div); return; }
    int n_eq = 0;
    for (int i = 0; i < n; i++) {
        if ((ps->particles[i].type == FLUID_PARTICLE && !ps->particles[i].on_surface) ||
             ps->particles[i].type == WALL_PARTICLE) {
            eq_idx[i] = n_eq++;
        } else {
            eq_idx[i] = -1;
        }
    }

    if (n_eq == 0) {
        for (int i = 0; i < n; i++) ps->particles[i].pressure = 0.0;
        free(eq_idx); free(vel_div);
        return;
    }

    /* ---- CSR行列の構築 ----
     *
     * row_nnz[] の計算は各行が独立 → 並列化可能。
     * row_ptr の prefix-sum は逐次。
     * 要素書き込みは各行が独立したメモリ領域に書くが、
     * fill[] カウンタへのアクセスが行ごとに独立しているため並列化可能。
     */
    int *row_nnz = (int *)calloc(n_eq, sizeof(int));
    if (!row_nnz) { free(eq_idx); free(vel_div); return; }

    /* row_nnz[] の計算 (並列) */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        if (eq_idx[i] < 0) continue;
        int ei = eq_idx[i];
        int cnt = 1; /* 対角 */
        for (int k = 0; k < nl->count[i]; k++) {
            int j = nl->neighbors[i * max_nb + k];
            if (ps->particles[j].type == DUMMY_PARTICLE) continue;
            if (eq_idx[j] >= 0) cnt++;
        }
        row_nnz[ei] = cnt;
    }

    /* prefix-sum で row_ptr を構築 (逐次) */
    int nnz = 0;
    int *row_ptr_tmp = (int *)malloc((n_eq + 1) * sizeof(int));
    if (!row_ptr_tmp) { free(row_nnz); free(eq_idx); free(vel_div); return; }
    row_ptr_tmp[0] = 0;
    for (int ei = 0; ei < n_eq; ei++) {
        nnz += row_nnz[ei];
        row_ptr_tmp[ei + 1] = nnz;
    }
    free(row_nnz);

    CsrMatrix *A = csr_alloc(n_eq, nnz);
    double    *b_vec = (double *)calloc(n_eq, sizeof(double));
    double    *x_vec = (double *)calloc(n_eq, sizeof(double));
    int       *fill  = (int    *)malloc(n_eq * sizeof(int));

    if (!A || !b_vec || !x_vec || !fill) {
        fprintf(stderr, "Error: pressure solver memory allocation failed "
                        "(n_eq=%d, nnz=%d)\n", n_eq, nnz);
        csr_free(A); free(b_vec); free(x_vec); free(fill);
        free(row_ptr_tmp); free(eq_idx); free(vel_div);
        return;
    }
    memcpy(A->row_ptr, row_ptr_tmp, (n_eq + 1) * sizeof(int));
    free(row_ptr_tmp);

    /* fill[] を row_ptr[ei] で初期化 (並列) */
    #pragma omp parallel for schedule(static)
    for (int ei = 0; ei < n_eq; ei++) fill[ei] = A->row_ptr[ei];

    /* 要素の書き込み (並列)
     * 各粒子 i は eq_idx[i] が示す行 ei に対してのみ書き込む。
     * fill[ei] は粒子ごとに独立しており、複数スレッドが同一の ei を
     * 書き込まないようにするため、eq_idx[i] が一意である性質を利用。
     * （eq_idx は単調増加の prefix-sum のため同一値は1回のみ出現）
     */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        if (eq_idx[i] < 0) continue;
        int ei = eq_idx[i];

        double sum_w = 0.0;
        for (int k = 0; k < nl->count[i]; k++) {
            int j = nl->neighbors[i * max_nb + k];
            if (ps->particles[j].type == DUMMY_PARTICLE) continue;
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                double diff = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += diff * diff;
            }
            double w = kernel_weight(sqrt(r2), re);
            sum_w += w;

            if (eq_idx[j] >= 0) {
                int pos = fill[ei]++;
                A->col_idx[pos] = eq_idx[j];
                A->val[pos]     = -coeff * w;
            }
        }

        /* 対角 */
        double diag_val = (g_config->ppe_type == 1)
                          ? g_config->c_ppe * coeff * sum_w
                          : coeff * sum_w;
        int pos = fill[ei]++;
        A->col_idx[pos] = ei;
        A->val[pos]     = diag_val;
        A->diag[ei]     = diag_val;

        /* 右辺 */
        if (ps->particles[i].type == WALL_PARTICLE)
            b_vec[ei] = 0.0;
        else if (g_config->hs_mode == 1)
            b_vec[ei] = source_term_hs(i, ps, nl, n0, dt, re);
        else
            b_vec[ei] = source_term_standard(i, ps, n0, dt, vel_div);
    }
    free(fill);

    /* 各行を列インデックスでソート (IC分解に必要)
     * 各行は独立したメモリ領域 → 並列化可能。
     */
    #pragma omp parallel for schedule(dynamic, 64)
    for (int ei = 0; ei < n_eq; ei++) {
        int lo = A->row_ptr[ei], hi = A->row_ptr[ei + 1];
        for (int p = lo + 1; p < hi; p++) {
            int    ci = A->col_idx[p];
            double vi = A->val[p];
            int q = p - 1;
            while (q >= lo && A->col_idx[q] > ci) {
                A->col_idx[q + 1] = A->col_idx[q];
                A->val[q + 1]     = A->val[q];
                q--;
            }
            A->col_idx[q + 1] = ci;
            A->val[q + 1]     = vi;
        }
    }

    /* ソルバー */
    int iters;
    if (g_config->solver_type == 1)
        iters = solve_iccg_sparse(A, x_vec, b_vec,
                                  g_config->cg_max_iter, g_config->cg_tolerance);
    else
        iters = solve_cg_sparse(A, x_vec, b_vec,
                                g_config->cg_max_iter, g_config->cg_tolerance);
    (void)iters;

    /* 結果を粒子に反映 (並列) */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++)
        ps->particles[i].pressure = (eq_idx[i] >= 0) ? x_vec[eq_idx[i]] : 0.0;

    csr_free(A);
    free(b_vec); free(x_vec); free(eq_idx); free(vel_div);
}
