#include <math.h>
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
 *   <∇P>_i = (d / n0) * Σ_{j≠i} [(P_j - P_min) / |r_j - r_i|^2]
 *            * (r_j - r_i) * w(|r_j - r_i|, re)
 *
 * P_min: 粒子i及びその近傍の最小圧力（引張不安定性対策）
 *
 * 補正加速度 = -(1/ρ) * <∇P> をaccに格納
 */
void calc_pressure_gradient(ParticleSystem *ps, NeighborList *nl)
{
    double re = g_config->influence_radius_lap;
    double n0 = ps->n0;
    double grad_coeff = (double)DIM / n0;

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        /* 近傍粒子（自分を含む）の最小圧力を求める（引張不安定性対策） */
        double p_min = ps->particles[i].pressure;
        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
            if (ps->particles[j].pressure < p_min)
                p_min = ps->particles[j].pressure;
        }

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
            if (r2 < 1.0e-20) continue;

            double r = sqrt(r2);
            double w = kernel_weight(r, re);
            double dp = ps->particles[j].pressure - p_min;

            for (int d = 0; d < DIM; d++) {
                grad[d] += dp / r2 * dr[d] * w;
            }
        }

        /* acc = -(1/ρ) * ∇P */
        for (int d = 0; d < DIM; d++) {
            ps->particles[i].acc[d] = -grad_coeff * grad[d] / g_config->density;
        }
    }
}

/*
 * 自由表面判定
 *   粒子数密度 n_i < threshold * n0 なら自由表面とみなす
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
