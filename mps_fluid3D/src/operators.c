#include <math.h>
#include <stdlib.h>
#include "operators.h"
#include "kernel.h"
#include "sim_config.h"

/*
 * 粒子数密度の計算
 *   n_i = Σ_{j≠i} w(|r_j - r_i|, re)
 */
void calc_particle_number_density(ParticleSystem *ps, NeighborList *nl)
{
    double re = g_config->influence_radius_n;

    for (int i = 0; i < ps->num; i++) {
        double ni = 0.0;
        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
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
 */
void calc_viscosity_term(ParticleSystem *ps, NeighborList *nl)
{
    double re = g_config->influence_radius_lap;
    double n0 = ps->n0;
    double lambda = ps->lambda;
    double coeff = 2.0 * DIM / (n0 * lambda);

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double lap[DIM];
        for (int d = 0; d < DIM; d++) lap[d] = 0.0;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
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
            ps->particles[i].acc[d] += g_config->viscosity * coeff * lap[d];
        }
    }
}

/*
 * 圧力勾配（勾配モデル）
 *
 * [標準モード: cmps_gradient=0]
 *   dp = P_j - P_min
 *   P_min: 粒子i及びその近傍の最小圧力（引張不安定性対策）
 *
 * [CMPSモード: cmps_gradient=1]
 *   dp = P_i + P_j - P_i_min - P_j_min  （対称型）
 *   P_i_min: 粒子iの近傍最小圧力
 *   P_j_min: 粒子jの近傍最小圧力
 *
 * 補正加速度 = -(1/ρ) * (d/n0) * Σ [dp/r^2 * dr * w] をaccに格納
 */
void calc_pressure_gradient(ParticleSystem *ps, NeighborList *nl)
{
    double re = g_config->influence_radius_lap;
    double n0 = ps->n0;
    double grad_coeff = (double)DIM / n0;
    double l0 = g_config->particle_distance;
    int cmps = g_config->cmps_gradient;

    /* 全粒子の近傍最小圧力を事前計算
     * 標準モードでは p_min_arr[i] のみ使用、CMPSモードでは p_min_arr[j] も使用 */
    double *p_min_arr = malloc(ps->num * sizeof(double));
    if (!p_min_arr) return;

    for (int i = 0; i < ps->num; i++) {
        double p_min = ps->particles[i].pressure;
        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
            if (ps->particles[j].pressure < p_min)
                p_min = ps->particles[j].pressure;
        }
        p_min_arr[i] = p_min;
    }

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double pi_min = p_min_arr[i];

        double grad[DIM];
        for (int d = 0; d < DIM; d++) grad[d] = 0.0;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
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
            if (cmps) {
                /* CMPSモード: 対称型 */
                dp = ps->particles[i].pressure + ps->particles[j].pressure
                     - pi_min - p_min_arr[j];
            } else {
                /* 標準モード */
                dp = ps->particles[j].pressure - pi_min;
            }

            for (int d = 0; d < DIM; d++) {
                grad[d] += dp / r2 * dr[d] * w;
            }
        }

        /* acc = -(1/ρ) * ∇P */
        for (int d = 0; d < DIM; d++) {
            ps->particles[i].acc[d] = -grad_coeff * grad[d] / g_config->density;
        }
    }

    free(p_min_arr);
}

/*
 * 負圧クランプ
 *   流体粒子の圧力を max(P_i, 0) に制限する。
 *   引張不安定性（粒子の凝集）を抑制する効果がある。
 *   params.txt の clamp_negative_pressure=1 のとき solve_pressure 直後に呼ぶ。
 */
void clamp_pressure(ParticleSystem *ps)
{
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == FLUID_PARTICLE &&
            ps->particles[i].pressure < 0.0) {
            ps->particles[i].pressure = 0.0;
        }
    }
}

/*
 * 自由表面判定 (粒子数密度ベース, 既存手法)
 *   n_i < threshold * n0 なら自由表面とみなす
 */
void judge_free_surface(ParticleSystem *ps, double threshold)
{
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == FLUID_PARTICLE) {
            ps->particles[i].on_surface = (ps->particles[i].n < threshold * ps->n0) ? 1 : 0;
        } else {
            ps->particles[i].on_surface = 0;
        }
    }
}

/*
 * 近傍粒子数の計算 (Natsui法用)
 *   influence_radius_n 内の粒子数を Ni として各粒子に記録する
 */
void calc_neighbor_count(ParticleSystem *ps, NeighborList *nl)
{
    double re2 = g_config->influence_radius_n * g_config->influence_radius_n;

    for (int i = 0; i < ps->num; i++) {
        int count = 0;
        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
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
 *   Ni < beta * N0 なら自由表面とみなす
 */
void judge_free_surface_by_count(ParticleSystem *ps, double beta)
{
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == FLUID_PARTICLE) {
            ps->particles[i].on_surface =
                (ps->particles[i].neighbor_count < beta * ps->n0_count) ? 1 : 0;
        } else {
            ps->particles[i].on_surface = 0;
        }
    }
}

/*
 * ポテンシャル型表面張力 (Colagrossi型)
 *   f(r) = -(r - l0) * (r - re_st)
 *   coeff = C_LL / (ρ * l0³)
 *   accに加算（陽的ステップ）
 */
void calc_surface_tension(ParticleSystem *ps, NeighborList *nl)
{
    double re_st = g_config->influence_radius_st;
    double l0    = g_config->particle_distance;
    double coeff = ps->C_LL / (g_config->density * l0 * l0 * l0);

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
            if (ps->particles[j].type != FLUID_PARTICLE) continue;

            double dr[DIM], r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                dr[d] = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += dr[d] * dr[d];
            }
            double r = sqrt(r2);
            if (r < 1.0e-9 * l0 || r >= re_st) continue;

            /* ポテンシャルの負の微分 = 粒子間力
             * r < l0 → 斥力, l0 < r < re_st → 引力 */
            double force_mag = -(r - l0) * (r - re_st);

            for (int d = 0; d < DIM; d++)
                ps->particles[i].acc[d] += coeff * force_mag * (dr[d] / r);
        }
    }
}

/*
 * 粒子間衝突モデル（越塚 2003 に基づく）
 *
 * 粒子間距離が collision_dist = collision_distance_ratio * l0 を下回り、かつ接近中
 * （相対速度の法線成分 > 0）の場合に衝突インパルスを適用する。
 *
 * mi = mj（同一密度）のとき、インパルス式は密度がキャンセルされ
 *   Δv_i = -(1+e)/2 * v_rel_n * n_ij
 * となる（n_ij: i→j の単位ベクトル）。
 *
 * 呼び出しタイミング: 陽的ステップ（仮位置更新）直後、近傍探索前。
 */
void collision(ParticleSystem *ps, NeighborList *nl)
{
    double l0        = g_config->particle_distance;
    double col_dist  = g_config->collision_distance_ratio * l0;
    double col_dist2 = col_dist * col_dist;
    double e         = g_config->restitution_coeff;
    double dt        = g_config->dt;

    /* 衝突後速度の一時バッファ（元の速度で初期化） */
    double (*vel_after)[DIM] = malloc(ps->num * sizeof(*vel_after));
    if (!vel_after) return;

    for (int i = 0; i < ps->num; i++) {
        for (int d = 0; d < DIM; d++)
            vel_after[i][d] = ps->particles[i].vel[d];
    }

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double vi[DIM];
        for (int d = 0; d < DIM; d++)
            vi[d] = ps->particles[i].vel[d];

        /* 近傍リストのみ走査 (col_dist < re なので衝突ペアは全て含まれる) */
        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
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

    /* 速度・位置の更新 */
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;
        for (int d = 0; d < DIM; d++) {
            double dv = vel_after[i][d] - ps->particles[i].vel[d];
            ps->particles[i].pos[d] += dt * dv;
            ps->particles[i].vel[d]  = vel_after[i][d];
        }
    }

    free(vel_after);
}
