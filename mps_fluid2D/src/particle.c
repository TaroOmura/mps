#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "particle.h"
#include "kernel.h"
#include "config.h"

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
 * 全流体粒子の粒子数密度を計算し、最大値を n0 とする
 * （最大値を持つ粒子が最も完全な近傍を持つ内部粒子）
 *
 *   n0 = max_i { Σ_{j≠i} w(|r_j - r_i|, re) }
 *   λ  = Σ_{j≠i} |r_j - r_i|^2 * w(|r_j - r_i|, re) / n0
 */
void particle_system_calc_initial_params(ParticleSystem *ps)
{
    double re = INFLUENCE_RADIUS;
    double max_n0 = 0.0;
    double max_lambda_num = 0.0;

    int found = 0;
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double n0_i = 0.0;
        double lambda_num_i = 0.0;

        for (int j = 0; j < ps->num; j++) {
            if (j == i) continue;
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                double diff = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += diff * diff;
            }
            double r = sqrt(r2);
            double w = kernel_weight(r, re);
            if (w > 0.0) {
                n0_i += w;
                lambda_num_i += r2 * w;
            }
        }

        if (n0_i > max_n0) {
            max_n0 = n0_i;
            max_lambda_num = lambda_num_i;
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
    ps->lambda = (max_n0 > 1.0e-10) ? max_lambda_num / max_n0 : 1.0;

    printf("Initial params: n0 = %.6f, lambda = %.6f\n", ps->n0, ps->lambda);
}
