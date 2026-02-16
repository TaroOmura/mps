#include <math.h>
#include "boundary.h"
#include "sim_config.h"

/* 壁粒子の速度・加速度をゼロに固定 */
void apply_wall_boundary(ParticleSystem *ps)
{
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == WALL_PARTICLE ||
            ps->particles[i].type == GHOST_PARTICLE) {
            for (int d = 0; d < DIM; d++) {
                ps->particles[i].vel[d] = 0.0;
                ps->particles[i].acc[d] = 0.0;
            }
        }
    }
}

/* 自由表面の圧力をゼロに設定 */
void apply_surface_pressure(ParticleSystem *ps)
{
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].on_surface) {
            ps->particles[i].pressure = 0.0;
        }
    }
}

/*
 * 粒子間衝突処理 (Koshizuka & Oka, 1996 に基づく)
 */
void handle_collision(ParticleSystem *ps, NeighborList *nl)
{
    double collision_dist = g_config->particle_distance * 0.5;
    double restitution = 0.2;

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);

            double dr[DIM];
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                dr[d] = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += dr[d] * dr[d];
            }
            double r = sqrt(r2);

            if (r >= collision_dist || r < 1.0e-10) continue;

            double vn = 0.0;
            for (int d = 0; d < DIM; d++) {
                vn += (ps->particles[j].vel[d] - ps->particles[i].vel[d]) * dr[d] / r;
            }

            if (vn >= 0.0) continue;

            double impulse = -(1.0 + restitution) * vn * 0.5;

            for (int d = 0; d < DIM; d++) {
                double dv = impulse * dr[d] / r;
                ps->particles[i].vel[d] -= dv;
                if (ps->particles[j].type == FLUID_PARTICLE) {
                    ps->particles[j].vel[d] += dv;
                }
            }
        }
    }
}

/*
 * 壁面反発力
 */
void apply_wall_repulsion(ParticleSystem *ps, NeighborList *nl)
{
    double l0 = g_config->particle_distance;
    double dt = g_config->dt;
    double coeff = g_config->wall_repulsion_coeff;

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        for (int k = 0; k < nl->count[i]; k++) {
            int j = neighbor_get(nl, i, k);
            if (ps->particles[j].type != WALL_PARTICLE) continue;

            double dr[DIM];
            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                dr[d] = ps->particles[i].pos[d] - ps->particles[j].pos[d];
                r2 += dr[d] * dr[d];
            }
            double r = sqrt(r2);

            if (r >= l0 || r < 1.0e-10) continue;

            double overlap = 1.0 - r / l0;
            double force_mag = coeff * overlap * overlap;

            for (int d = 0; d < DIM; d++) {
                ps->particles[i].vel[d] += dt * force_mag * dr[d] / r;
            }
        }
    }
}

/*
 * 壁面位置クランプ (3D)
 * 上面は開放境界のためクランプしない。
 */
void clamp_to_walls(ParticleSystem *ps)
{
    double half_l0 = g_config->particle_distance * 0.5;
    double x_min = g_config->domain_min[0] + half_l0;
    double x_max = g_config->domain_max[0] - half_l0;
    double y_min = g_config->domain_min[1] + half_l0;
    double z_min = g_config->domain_min[2] + half_l0;
    double z_max = g_config->domain_max[2] - half_l0;
    double rest  = g_config->wall_restitution;

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        /* 左壁 */
        if (ps->particles[i].pos[0] < x_min) {
            ps->particles[i].pos[0] = x_min;
            if (ps->particles[i].vel[0] < 0.0)
                ps->particles[i].vel[0] *= -rest;
        }
        /* 右壁 */
        if (ps->particles[i].pos[0] > x_max) {
            ps->particles[i].pos[0] = x_max;
            if (ps->particles[i].vel[0] > 0.0)
                ps->particles[i].vel[0] *= -rest;
        }
        /* 底面 */
        if (ps->particles[i].pos[1] < y_min) {
            ps->particles[i].pos[1] = y_min;
            if (ps->particles[i].vel[1] < 0.0)
                ps->particles[i].vel[1] *= -rest;
        }
        /* 前壁 */
        if (ps->particles[i].pos[2] < z_min) {
            ps->particles[i].pos[2] = z_min;
            if (ps->particles[i].vel[2] < 0.0)
                ps->particles[i].vel[2] *= -rest;
        }
        /* 後壁 */
        if (ps->particles[i].pos[2] > z_max) {
            ps->particles[i].pos[2] = z_max;
            if (ps->particles[i].vel[2] > 0.0)
                ps->particles[i].vel[2] *= -rest;
        }
    }
}

/* 計算領域外の粒子をゴースト粒子に変更 */
void remove_out_of_bounds(ParticleSystem *ps, double xmin, double xmax,
                          double ymin, double ymax,
                          double zmin, double zmax)
{
    double margin = g_config->particle_distance * (g_config->wall_layers + 1);
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double x = ps->particles[i].pos[0];
        double y = ps->particles[i].pos[1];
        double z = ps->particles[i].pos[2];

        if (x < xmin - margin || x > xmax + margin ||
            y < ymin - margin || y > ymax + margin ||
            z < zmin - margin || z > zmax + margin ||
            isnan(x) || isnan(y) || isnan(z)) {
            ps->particles[i].type = GHOST_PARTICLE;
            for (int d = 0; d < DIM; d++) {
                ps->particles[i].vel[d] = 0.0;
                ps->particles[i].acc[d] = 0.0;
            }
        }
    }
}
