/* OpenMP parallelization added */
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "operators.h"
#include "kernel.h"
#include "sim_config.h"

/*
 * 粒子数密度の計算
 *   n_i = Σ_{j≠i} w(|r_j - r_i|, re)
 *
 * OpenMP並列化:
 *   各粒子iは独立して n_i を計算し ps->particles[i].n に書き込む。
 *   近傍リストは読み取り専用 → データ競合なし。
 */
void calc_particle_number_density(ParticleSystem *ps, NeighborList *nl)
{
    double re = g_config->influence_radius_n;
    int n = ps->num;
    int max_nb = nl->max_neighbors;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        double ni = 0.0;
        for (int k = 0; k < nl->count[i]; k++) {
            int j = nl->neighbors[i * max_nb + k];
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                double diff = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += diff * diff;
            }
            double r = sqrt(r2);
            ni += kernel_weight(r, re);
        }
        ps->particles[i].n = ni;
    }
}

/*
 * 粘性項（ラプラシアンモデル）
 *   <∇²u>_i = (2d / (n0 * λ)) * Σ_{j≠i} (u_j - u_i) * w(|r_j - r_i|, re)
 *
 * 粘性力 = ν * <∇²u> を加速度に加算
 *
 * OpenMP並列化:
 *   各粒子iの acc[d] への加算は粒子iのみが書き込む。
 *   近傍粒子jの vel[] は読み取り専用 → データ競合なし。
 */
void calc_viscosity_term(ParticleSystem *ps, NeighborList *nl)
{
    double re = g_config->influence_radius_lap;
    double n0 = ps->n0;
    double lambda = ps->lambda;
    double coeff = 2.0 * DIM / (n0 * lambda);
    double visc = g_config->viscosity;
    int n = ps->num;
    int max_nb = nl->max_neighbors;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double lap[DIM];
        for (int d = 0; d < DIM; d++) lap[d] = 0.0;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = nl->neighbors[i * max_nb + k];
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                double diff = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += diff * diff;
            }
            double r = sqrt(r2);
            double w = kernel_weight(r, re);
            for (int d = 0; d < DIM; d++) {
                lap[d] += (ps->particles[j].vel[d] - ps->particles[i].vel[d]) * w;
            }
        }

        for (int d = 0; d < DIM; d++) {
            ps->particles[i].acc[d] += visc * coeff * lap[d];
        }
    }
}

/*
 * 圧力勾配（勾配モデル）
 *
 * [標準モード: cmps_gradient=0]
 *   dp = P_j - P_min
 * [CMPSモード: cmps_gradient=1] (Khayyer & Gotoh, 2008)
 *   dp = P_i + P_j - P_i_min - P_j_min
 * [Oochiモード: cmps_gradient=2] (Oochi, 2010)
 *   dp = P_i + P_j  （対称型、P_min補正なし）
 *
 * OpenMP並列化:
 *   ステップ1: p_min_arr[] の計算 - 各粒子iは独立して書き込む。
 *   ステップ2: 勾配計算・acc[] への書き込み - 各粒子iが acc[i] に書き込むのみ。
 *   両ステップともデータ競合なし。
 */
void calc_pressure_gradient(ParticleSystem *ps, NeighborList *nl)
{
    double re = g_config->influence_radius_lap;
    double n0 = ps->n0;
    double grad_coeff = (double)DIM / n0;
    double l0 = g_config->particle_distance;
    double density = g_config->density;
    int cmps = g_config->cmps_gradient;
    int n = ps->num;
    int max_nb = nl->max_neighbors;

    /* Oochiモード (cmps==2) ではP_min不使用のためスキップ */
    double *p_min_arr = NULL;
    if (cmps != 2) {
        p_min_arr = malloc(n * sizeof(double));
        if (!p_min_arr) return;

        /* ステップ1: 各粒子の近傍最小圧力を計算 (並列) */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) {
            double p_min = ps->particles[i].pressure;
            for (int k = 0; k < nl->count[i]; k++) {
                int j = nl->neighbors[i * max_nb + k];
                if (ps->particles[j].pressure < p_min)
                    p_min = ps->particles[j].pressure;
            }
            p_min_arr[i] = p_min;
        }
    }

    /* ステップ2: 圧力勾配の計算 (並列)
     * p_min_arr[] は読み取り専用 (ステップ1完了後)、
     * acc[i] への書き込みはスレッド間で衝突しない。
     */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double pi_min = p_min_arr ? p_min_arr[i] : 0.0;

        double grad[DIM];
        for (int d = 0; d < DIM; d++) grad[d] = 0.0;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = nl->neighbors[i * max_nb + k];
            double r2 = 0.0;
            double dr[DIM];
            for (int d = 0; d < DIM; d++) {
                dr[d] = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += dr[d] * dr[d];
            }
            if (r2 < l0 * l0 * 1.0e-12) continue;

            double r = sqrt(r2);
            double w = kernel_weight(r, re);
            double dp;
            if (cmps == 1) {
                /* CMPSモード (Khayyer & Gotoh, 2008): 対称型 */
                dp = ps->particles[i].pressure + ps->particles[j].pressure
                     - pi_min - p_min_arr[j];
            } else if (cmps == 2) {
                /* Oochiモード (Oochi, 2010): 対称型、P_min補正なし */
                dp = ps->particles[i].pressure + ps->particles[j].pressure;
            } else {
                /* 標準モード */
                dp = ps->particles[j].pressure - pi_min;
            }

            for (int d = 0; d < DIM; d++) {
                grad[d] += dp / r2 * dr[d] * w;
            }
        }

        for (int d = 0; d < DIM; d++) {
            ps->particles[i].acc[d] = -grad_coeff * grad[d] / density;
        }
    }

    free(p_min_arr);
}

/*
 * 負圧クランプ
 *
 * OpenMP並列化: 各粒子は独立 → データ競合なし。
 */
void clamp_pressure(ParticleSystem *ps)
{
    int n = ps->num;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type == FLUID_PARTICLE &&
            ps->particles[i].pressure < 0.0) {
            ps->particles[i].pressure = 0.0;
        }
    }
}

/*
 * 自由表面判定 (粒子数密度ベース)
 *
 * OpenMP並列化: 各粒子は独立 → データ競合なし。
 */
void judge_free_surface(ParticleSystem *ps, double threshold)
{
    int n = ps->num;
    double n0 = ps->n0;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type == FLUID_PARTICLE) {
            ps->particles[i].on_surface = (ps->particles[i].n < threshold * n0) ? 1 : 0;
        } else {
            ps->particles[i].on_surface = 0;
        }
    }
}

/*
 * 近傍粒子数の計算 (Natsui法用)
 *
 * OpenMP並列化: 各粒子iは独立して neighbor_count に書き込む。
 */
void calc_neighbor_count(ParticleSystem *ps, NeighborList *nl)
{
    double re2 = g_config->influence_radius_n * g_config->influence_radius_n;
    int n = ps->num;
    int max_nb = nl->max_neighbors;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        int count = 0;
        for (int k = 0; k < nl->count[i]; k++) {
            int j = nl->neighbors[i * max_nb + k];
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                double diff = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += diff * diff;
            }
            if (r2 < re2) count++;
        }
        ps->particles[i].neighbor_count = count;
    }
}

/*
 * 自由表面判定 (近傍粒子数ベース, Natsui法)
 *
 * OpenMP並列化: 各粒子は独立 → データ競合なし。
 */
void judge_free_surface_by_count(ParticleSystem *ps, double beta)
{
    int n = ps->num;
    int n0_count = ps->n0_count;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type == FLUID_PARTICLE) {
            ps->particles[i].on_surface =
                (ps->particles[i].neighbor_count < beta * n0_count) ? 1 : 0;
        } else {
            ps->particles[i].on_surface = 0;
        }
    }
}

/*
 * ポテンシャル型表面張力 (Colagrossi型)
 *   f(r) = -(r - l0) * (r - re_st)
 *   coeff = C_LL / (ρ * l0³)
 *
 * OpenMP並列化:
 *   各粒子iの acc[i] への書き込みはスレッド間で衝突しない。
 *   近傍粒子jの pos[] は読み取り専用 → データ競合なし。
 *
 * 注意: 元のコードは acc[i] += ... のみで acc[j] は更新しない（非対称）。
 *       そのため各iのループが完全に独立しており、並列化は安全。
 */
void calc_surface_tension(ParticleSystem *ps, NeighborList *nl)
{
    double re_st = g_config->influence_radius_st;
    double l0    = g_config->particle_distance;
    double coeff = ps->C_LL / (g_config->density * l0 * l0 * l0);
    int n = ps->num;
    int max_nb = nl->max_neighbors;

    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = nl->neighbors[i * max_nb + k];
            if (ps->particles[j].type != FLUID_PARTICLE) continue;

            double dr[DIM], r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                dr[d] = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += dr[d] * dr[d];
            }
            double r = sqrt(r2);
            if (r < 1.0e-9 * l0 || r >= re_st) continue;

            double force_mag = -(r - l0) * (r - re_st);

            for (int d = 0; d < DIM; d++)
                ps->particles[i].acc[d] += coeff * force_mag * (dr[d] / r);
        }
    }
}

/*
 * 固液表面張力（濡れモデル）
 *   流体粒子 i と壁粒子 j のペアに対して液液と同じポテンシャルを C_SL で適用する。
 *   壁粒子は固定なので流体粒子の加速度のみ更新する。
 *   C_SL = 0.5*(1+cosθ)*C_LL  (Young-Dupré式)
 *
 * OpenMP並列化:
 *   acc[i] への書き込みのみ、データ競合なし。
 */
void calc_surface_tension_SL(ParticleSystem *ps, NeighborList *nl)
{
    if (ps->C_SL == 0.0) return;

    double re_st   = g_config->influence_radius_st;
    double l0      = g_config->particle_distance;
    double coeff   = ps->C_SL / (g_config->density * l0 * l0 * l0);
    int n          = ps->num;
    int max_nb     = nl->max_neighbors;

    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = nl->neighbors[i * max_nb + k];
            if (ps->particles[j].type != WALL_PARTICLE) continue;

            double dr[DIM], r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                dr[d] = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += dr[d] * dr[d];
            }
            double r = sqrt(r2);
            if (r < 1.0e-9 * l0 || r >= re_st) continue;

            double force_mag = -(r - l0) * (r - re_st);

            for (int d = 0; d < DIM; d++)
                ps->particles[i].acc[d] += coeff * force_mag * (dr[d] / r);
        }
    }
}

/*
 * 粒子間衝突モデル（越塚 2003）
 *
 * OpenMP並列化の方針:
 *   衝突モデルは「衝突後速度バッファ vel_after[]」を介して各粒子iの速度変化を
 *   個別に計算し、その後まとめて更新する設計。
 *
 *   ループ内では vel_after[i] に書き込み、ps->particles[j].vel[] を読み取るのみ。
 *   vel_after[j] への書き込みは行わない（元のコードの設計に従う）。
 *   → データ競合なし、並列化可能。
 *
 *   注意: この実装は「各粒子iが独立して衝突を処理する」近似であり、
 *         同一ペア(i,j)に対して両方向の処理が非同期に起こりうるが、
 *         元の逐次コードでも同様の近似（vel_after[j]は更新しない）を採用している。
 */
void collision(ParticleSystem *ps, NeighborList *nl)
{
    double l0        = g_config->particle_distance;
    double col_dist  = g_config->collision_distance_ratio * l0;
    double col_dist2 = col_dist * col_dist;
    double e         = g_config->restitution_coeff;
    double dt        = g_config->dt;
    int n = ps->num;
    int max_nb = nl->max_neighbors;

    double (*vel_after)[DIM] = malloc(n * sizeof(*vel_after));
    if (!vel_after) return;

    /* 初期化 (並列) */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        for (int d = 0; d < DIM; d++)
            vel_after[i][d] = ps->particles[i].vel[d];
    }

    /* 衝突処理 (並列)
     * vel_after[i] への書き込みのみ。ps->particles[j].vel[] は読み取り専用。
     * vel_after[j] は更新しない（元コードと同じ設計）。
     */
    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double vi[DIM];
        for (int d = 0; d < DIM; d++)
            vi[d] = ps->particles[i].vel[d];

        for (int k = 0; k < nl->count[i]; k++) {
            int j = nl->neighbors[i * max_nb + k];
            if (ps->particles[j].type == GHOST_PARTICLE) continue;

            double dr[DIM];
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                dr[d] = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += dr[d] * dr[d];
            }

            if (r2 >= col_dist2) continue;

            double r = sqrt(r2);
            if (r < 1.0e-12 * l0) continue;

            double rel_v_n = 0.0;
            for (int d = 0; d < DIM; d++)
                rel_v_n += (vi[d] - ps->particles[j].vel[d]) * (dr[d] / r);

            if (rel_v_n <= 0.0) continue;

            double impulse = (1.0 + e) * 0.5 * rel_v_n;
            for (int d = 0; d < DIM; d++)
                vi[d] -= impulse * (dr[d] / r);
        }

        for (int d = 0; d < DIM; d++)
            vel_after[i][d] = vi[d];
    }

    /* 速度・位置の更新 (並列) */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;
        for (int d = 0; d < DIM; d++) {
            double dv = vel_after[i][d] - ps->particles[i].vel[d];
            ps->particles[i].pos[d] += dt * dv;
            ps->particles[i].vel[d]  = vel_after[i][d];
        }
    }

    free(vel_after);
}
