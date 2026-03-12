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
 * [CMPSモード: cmps_gradient=1] (Khayyer & Gotoh, 2008)
 *   dp = P_i + P_j - P_i_min - P_j_min  （対称型）
 *   P_i_min: 粒子iの近傍最小圧力
 *   P_j_min: 粒子jの近傍最小圧力
 *
 * [Oochiモード: cmps_gradient=2] (Oochi, 2010)
 *   dp = P_i + P_j  （対称型、P_min補正なし）
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
     * 標準モードでは p_min_arr[i] のみ使用、CMPSモードでは p_min_arr[j] も使用
     * Oochiモード (cmps==2) ではP_min不使用のためスキップ */
    double *p_min_arr = NULL;
    if (cmps != 2) {
        p_min_arr = malloc(ps->num * sizeof(double));
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
    }

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double pi_min = p_min_arr ? p_min_arr[i] : 0.0;

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
 * 固液表面張力（濡れモデル）
 *   流体粒子 i と壁粒子 j のペアに対して液液と同じポテンシャルを C_SL で適用する。
 *   壁粒子は固定なので流体粒子の加速度のみ更新する。
 *   C_SL = 0.5*(1+cosθ)*C_LL  (Young-Dupré式)
 */
void calc_surface_tension_SL(ParticleSystem *ps, NeighborList *nl)
{
    if (ps->C_SL == 0.0) return;

    double re_st = g_config->influence_radius_st;
    double l0    = g_config->particle_distance;
    double coeff = ps->C_SL / (g_config->density * l0 * l0 * l0);

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
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
 * ハイブリッド摩擦モデル (Hattori & Koshizuka 2019, Mechanical Engineering Journal)
 *
 * 接触線粒子の判定:
 *   - WALL粒子が re_st 内に存在する
 *   - かつ流体粒子のみの局所粒子数密度 n_f < surface_threshold * n0
 *     (標準の on_surface は壁粒子込みの n を使うため接触線で誤検出するため使わない)
 *
 * 壁法線の推定 (平面・球面・任意曲面に対応):
 *   1. 流体粒子 i の近傍から最近傍壁粒子 j* を特定する
 *   2. j* の近傍リスト内の WALL粒子群を用いて j* における壁外向き法線を推定:
 *        n_w = Σ_{k=WALL in nl[j*]} w(r) * (j* - k)/|j* - k|  → 正規化
 *      → 平面壁: 隣接壁粒子が平面状に並ぶため垂直成分のみ残る
 *      → 球面壁: 放射方向に向く（任意曲面に対応）
 *
 * 自由表面外向き法線の推定:
 *   n_lv = -Σ_{j=FLUID in nl[i]} w(r) * (j - i)/|j - i|  → 正規化
 *   (液体内部方向の逆 = 気体側への外向き法線)
 *
 * 接触角の計算:
 *   cos(θ_i) = -n_w · n_lv
 *
 * 滑り方向ベクトル (n_lv の壁接線面への射影):
 *   e_slide = (n_lv - (n_lv · n_w) * n_w)  → 正規化
 *
 * ハイブリッド切り替え (v_c = friction_velocity_threshold):
 *   |v_i| <= v_c → 静摩擦: |θ_i - θ_s| <= Δθ_c なら vel[i] = 0
 *   |v_i| >  v_c → 動摩擦: 滑り速度成分を直接減衰させる (vel を直接修正)
 *                    Δv_n = -dt · α · ν · (2D/n0λ) · Σ_{j=WALL} w(r_ij) · v_n
 *                    v_n = v_i · e_slide (滑り速度成分)
 *                    注: acc は implicit ステップで上書きされるため vel に直接作用する
 */
void apply_friction(ParticleSystem *ps, NeighborList *nl)
{
    if (!g_config->friction_enabled) return;

    double re_st    = g_config->influence_radius_st;
    double l0       = g_config->particle_distance;
    double theta_s  = g_config->wetting_angle_SL    * M_PI / 180.0;
    double delta_th = g_config->friction_delta_theta * M_PI / 180.0;
    double v_c      = g_config->friction_velocity_threshold;
    double alpha    = g_config->dynamic_friction_alpha;
    double nu       = g_config->viscosity;
    double dyn_coeff = alpha * nu * 2.0 * DIM / (ps->n0 * ps->lambda);

    double re_n      = g_config->influence_radius_n;
    double surf_thr  = g_config->surface_threshold;

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        /* --- Step 1: 接触線検出 ---
         * 標準の on_surface フラグは壁粒子込みの n で判定するため、
         * 液滴底部（壁粒子が下にある）では n ≈ n0 となり on_surface=0 になる。
         * ここでは (a) 流体粒子のみの局所粒子数密度 n_f が surf_thr*n0 未満
         *          (b) かつ WALL粒子が re_st 内に存在
         * の両方を満たす粒子を接触線粒子とする。 */
        double n_fluid = 0.0;
        int    j_star  = -1;
        double r2_star = 1.0e30;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                double diff = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += diff * diff;
            }
            if (ps->particles[j].type == FLUID_PARTICLE) {
                if (r2 < re_n * re_n)
                    n_fluid += kernel_weight(sqrt(r2), re_n);
            } else if (ps->particles[j].type == WALL_PARTICLE) {
                if (r2 < re_st * re_st && r2 < r2_star) {
                    r2_star = r2; j_star = j;
                }
            }
        }
        if (j_star < 0) continue;                        /* 壁近傍でない */
        if (n_fluid >= surf_thr * ps->n0) continue;      /* 流体内部粒子 */

        /* --- Step 2: j* の近傍 WALL粒子群から壁法線 n_w を推定 --- */
        double n_w[DIM] = {0.0, 0.0, 0.0};
        for (int k = 0; k < nl->count[j_star]; k++) {
            int m = neighbor_get(nl, j_star, k);
            if (ps->particles[m].type != WALL_PARTICLE) continue;
            double dr[DIM], r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                dr[d] = ps->particles[j_star].pos[d] - ps->particles[m].pos[d];
                r2 += dr[d] * dr[d];
            }
            double r = sqrt(r2);
            if (r < 1.0e-9 * l0) continue;
            double w = kernel_weight(r, re_st);
            for (int d = 0; d < DIM; d++)
                n_w[d] += w * dr[d] / r;
        }
        double n_w_mag = 0.0;
        for (int d = 0; d < DIM; d++) n_w_mag += n_w[d] * n_w[d];
        n_w_mag = sqrt(n_w_mag);
        if (n_w_mag < 1.0e-9) continue;
        for (int d = 0; d < DIM; d++) n_w[d] /= n_w_mag;

        /* --- Step 3: nl[i] の FLUID粒子から自由表面外向き法線 n_lv を推定 --- */
        double n_lv[DIM] = {0.0, 0.0, 0.0};
        int fluid_count = 0;
        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
            if (ps->particles[j].type != FLUID_PARTICLE) continue;
            double dr[DIM], r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                dr[d] = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += dr[d] * dr[d];
            }
            double r = sqrt(r2);
            if (r < 1.0e-9 * l0) continue;
            double w = kernel_weight(r, re_st);
            for (int d = 0; d < DIM; d++)
                n_lv[d] -= w * dr[d] / r;  /* 液体内部方向の逆 = 外向き */
            fluid_count++;
        }
        if (fluid_count == 0) continue;
        double n_lv_mag = 0.0;
        for (int d = 0; d < DIM; d++) n_lv_mag += n_lv[d] * n_lv[d];
        n_lv_mag = sqrt(n_lv_mag);
        if (n_lv_mag < 1.0e-9) continue;
        for (int d = 0; d < DIM; d++) n_lv[d] /= n_lv_mag;

        /* --- Step 4: 接触角 θ_i --- */
        double cos_th = 0.0;
        for (int d = 0; d < DIM; d++)
            cos_th -= n_w[d] * n_lv[d];  /* cos(θ) = -n_w · n_lv */
        cos_th = fmax(-1.0, fmin(1.0, cos_th));
        double theta_i = acos(cos_th);

        /* --- Step 5: 滑り方向ベクトル e_slide (n_lv の壁接線面射影) --- */
        double dot_lv_w = 0.0;
        for (int d = 0; d < DIM; d++)
            dot_lv_w += n_lv[d] * n_w[d];
        double e_slide[DIM];
        double e_slide_mag = 0.0;
        for (int d = 0; d < DIM; d++) {
            e_slide[d] = n_lv[d] - dot_lv_w * n_w[d];
            e_slide_mag += e_slide[d] * e_slide[d];
        }
        e_slide_mag = sqrt(e_slide_mag);
        if (e_slide_mag < 1.0e-9) continue;
        for (int d = 0; d < DIM; d++) e_slide[d] /= e_slide_mag;

        /* --- Step 6: 速度の大きさで静/動を切り替え --- */
        double v2 = 0.0;
        for (int d = 0; d < DIM; d++)
            v2 += ps->particles[i].vel[d] * ps->particles[i].vel[d];

        if (sqrt(v2) <= v_c) {
            /* 静摩擦: 接触角偏差が閾値以内なら仮速度をゼロに固定 */
            if (fabs(theta_i - theta_s) <= delta_th) {
                for (int d = 0; d < DIM; d++)
                    ps->particles[i].vel[d] = 0.0;
            }
        } else {
            /* 動摩擦: 滑り方向速度成分を直接減衰させる
             *   Δv_n = -dt · α · ν · (2D/n0λ) · Σ_{j=WALL} w(r_ij) · v_n
             * acc は implicit ステップで上書きされるため vel を直接修正する */
            double sum_w_wall = 0.0;
            for (int k = 0; k < nl->count[i]; k++) {
                int j = neighbor_get(nl, i, k);
                if (ps->particles[j].type != WALL_PARTICLE) continue;
                double r2 = 0.0;
                for (int d = 0; d < DIM; d++) {
                    double diff = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                    r2 += diff * diff;
                }
                double r = sqrt(r2);
                if (r < 1.0e-9 * l0) continue;
                sum_w_wall += kernel_weight(r, re_st);
            }
            double v_n = 0.0;
            for (int d = 0; d < DIM; d++)
                v_n += ps->particles[i].vel[d] * e_slide[d];
            /* Δv_n: 符号をv_nと逆に、大きさはmin(|dyn_coeff*dt*Σw*v_n|, |v_n|) */
            double dv_n = -dyn_coeff * g_config->dt * sum_w_wall * v_n;
            if (dv_n >  v_n) dv_n =  v_n;   /* 逆転防止 */
            if (dv_n < -v_n) dv_n = -v_n;
            for (int d = 0; d < DIM; d++)
                ps->particles[i].vel[d] += dv_n * e_slide[d];
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
