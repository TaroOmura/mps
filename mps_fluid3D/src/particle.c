#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "particle.h"
#include "kernel.h"
#include "sim_config.h"

ParticleSystem *particle_system_create(int capacity)
{
    ParticleSystem *ps = (ParticleSystem *)malloc(sizeof(ParticleSystem));
    if (!ps) return NULL;

    ps->particles = (Particle *)calloc(capacity, sizeof(Particle));
    if (!ps->particles) {
        free(ps);
        return NULL;
    }
    ps->num = 0;
    ps->capacity = capacity;
    ps->n0 = 0.0;
    ps->lambda = 0.0;
    return ps;
}

void particle_system_free(ParticleSystem *ps)
{
    if (!ps) return;
    free(ps->particles);
    free(ps);
}

int particle_system_add(ParticleSystem *ps, double pos[], double vel[], int type)
{
    if (ps->num >= ps->capacity) {
        fprintf(stderr, "Error: particle capacity exceeded (%d)\n", ps->capacity);
        return -1;
    }

    Particle *p = &ps->particles[ps->num];
    for (int d = 0; d < DIM; d++) {
        p->pos[d] = pos[d];
        p->vel[d] = vel[d];
        p->acc[d] = 0.0;
    }
    p->pressure = 0.0;
    p->n = 0.0;
    p->type = type;
    p->on_surface = 0;

    ps->num++;
    return ps->num - 1;
}

/*
 * 初期粒子配置から基準粒子数密度 n0 とラプラシアンモデル用パラメータ λ を計算
 *
 * n0     : calc_particle_number_density と同じ influence_radius_n で計算
 * lambda : ラプラシアンモデルの固有値補正係数。influence_radius_lap で計算
 *
 * 両者に異なる半径を使う場合、n0 と n_i は必ず同じ半径で計算しなければ
 * 自由表面判定および PPE 右辺 (n_i - n0)/n0 が破綻する。
 */
void particle_system_calc_initial_params(ParticleSystem *ps)
{
    double re_n   = g_config->influence_radius_n;
    double re_lap = g_config->influence_radius_lap;

    double max_n0        = 0.0;
    double max_lambda_num = 0.0;
    double max_lambda_den = 0.0;

    int found = 0;
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double n0_i        = 0.0;
        double lambda_num_i = 0.0;
        double lambda_den_i = 0.0;

        for (int j = 0; j < ps->num; j++) {
            if (j == i) continue;
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                double diff = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += diff * diff;
            }
            double r = sqrt(r2);

            /* n0: influence_radius_n で計算 (calc_particle_number_density と同じ) */
            double wn = kernel_weight(r, re_n);
            if (wn > 0.0)
                n0_i += wn;

            /* lambda: influence_radius_lap で計算 */
            double wl = kernel_weight(r, re_lap);
            if (wl > 0.0) {
                lambda_num_i += r2 * wl;
                lambda_den_i += wl;
            }
        }

        if (n0_i > max_n0) {
            max_n0        = n0_i;
            max_lambda_num = lambda_num_i;
            max_lambda_den = lambda_den_i;
            found = 1;
        }
    }

    if (!found) {
        fprintf(stderr, "Error: no fluid particle found for initial params\n");
        ps->n0 = 1.0;
        ps->lambda = 1.0;
        return;
    }

    ps->n0 = max_n0;
    if (g_config->use_analytical_lambda) {
        ps->lambda = re_lap * re_lap * (double)DIM * (double)(DIM - 1)
                     / ((double)(DIM + 1) * (double)(DIM + 2));
    } else {
        ps->lambda = (max_lambda_den > 1.0e-10) ? max_lambda_num / max_lambda_den : 1.0;
    }

    printf("Initial params: n0 = %.6f (re_n=%.4f)  lambda = %.6f (re_lap=%.4f)%s\n",
           ps->n0, re_n, ps->lambda, re_lap,
           g_config->use_analytical_lambda ? "  [analytical]" : "");
}
