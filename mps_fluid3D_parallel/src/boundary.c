/* OpenMP parallelization added */
#include <math.h>
#include <omp.h>
#include "boundary.h"
#include "sim_config.h"

/*
 * 壁粒子の速度・加速度をゼロに固定
 *
 * OpenMP並列化: 各粒子は独立して書き込む → データ競合なし。
 */
void apply_wall_boundary(ParticleSystem *ps)
{
    int n = ps->num;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type == WALL_PARTICLE ||
            ps->particles[i].type == GHOST_PARTICLE) {
            for (int d = 0; d < DIM; d++) {
                ps->particles[i].vel[d] = 0.0;
                ps->particles[i].acc[d] = 0.0;
            }
        }
    }
}

/*
 * 計算領域外の粒子をゴースト粒子に変更
 *
 * OpenMP並列化: 各粒子は独立して書き込む → データ競合なし。
 */
void remove_out_of_bounds(ParticleSystem *ps, double xmin, double xmax,
                          double ymin, double ymax,
                          double zmin, double zmax)
{
    double margin = g_config->particle_distance * (g_config->wall_layers + 1);
    int n = ps->num;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double x = ps->particles[i].pos[0];
        double y = ps->particles[i].pos[1];
        double z = ps->particles[i].pos[2];

        if (x < xmin - margin || x > xmax + margin ||
            y < ymin - margin || y > ymax + margin ||
            z < zmin - margin || z > zmax + margin ||
            isnan(x) || isnan(y) || isnan(z)) {
            ps->particles[i].type = GHOST_PARTICLE;
            ps->particles[i].pressure = 0.0;
            ps->particles[i].n = 0.0;
            for (int d = 0; d < DIM; d++) {
                ps->particles[i].vel[d] = 0.0;
                ps->particles[i].acc[d] = 0.0;
            }
        }
    }
}
