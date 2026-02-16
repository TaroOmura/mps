#include <math.h>
#include "boundary.h"
#include "config.h"

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
 *
 * 2粒子が衝突距離以下に接近し、かつ互いに近づいている場合、
 * 接近方向の相対速度を反転・減衰させる。
 *
 *   衝突距離: l0 * 0.5
 *   反発係数: 0.2 (非弾性衝突)
 */
void handle_collision(ParticleSystem *ps, NeighborList *nl)
{
    double collision_dist = PARTICLE_DISTANCE * 0.5;
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

            /* 接線方向の単位ベクトル n = dr / r */
            /* 相対速度の法線成分: vn = (v_j - v_i) · n */
            double vn = 0.0;
            for (int d = 0; d < DIM; d++) {
                vn += (ps->particles[j].vel[d] - ps->particles[i].vel[d]) * dr[d] / r;
            }

            /* 接近中 (vn < 0) の場合のみ衝突処理 */
            if (vn >= 0.0) continue;

            /* 速度補正量 */
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
 * 流体粒子 i と壁粒子 j の距離 r < l0 のとき、
 * overlap に比例した反発加速度を流体粒子に与える。
 */
void apply_wall_repulsion(ParticleSystem *ps, NeighborList *nl)
{
    double l0 = PARTICLE_DISTANCE;

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
            double force_mag = WALL_REPULSION_COEFF * overlap * overlap;

            /* 反発加速度を直接速度に反映 */
            for (int d = 0; d < DIM; d++) {
                ps->particles[i].vel[d] += DT * force_mag * dr[d] / r;
            }
        }
    }
}

/*
 * 壁面位置クランプ
 * 壁ジオメトリを直接使い、貫通した粒子を押し戻す最後の砦。
 * 上面は開放境界のためクランプしない。
 */
void clamp_to_walls(ParticleSystem *ps)
{
    double half_l0 = PARTICLE_DISTANCE * 0.5;
    double x_min = DOMAIN_X_MIN + half_l0;
    double x_max = DOMAIN_X_MAX - half_l0;
    double y_min = DOMAIN_Y_MIN + half_l0;

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        /* 左壁 */
        if (ps->particles[i].pos[0] < x_min) {
            ps->particles[i].pos[0] = x_min;
            if (ps->particles[i].vel[0] < 0.0)
                ps->particles[i].vel[0] *= -WALL_RESTITUTION;
        }
        /* 右壁 */
        if (ps->particles[i].pos[0] > x_max) {
            ps->particles[i].pos[0] = x_max;
            if (ps->particles[i].vel[0] > 0.0)
                ps->particles[i].vel[0] *= -WALL_RESTITUTION;
        }
        /* 底面 */
        if (ps->particles[i].pos[1] < y_min) {
            ps->particles[i].pos[1] = y_min;
            if (ps->particles[i].vel[1] < 0.0)
                ps->particles[i].vel[1] *= -WALL_RESTITUTION;
        }
    }
}

/* 計算領域外の粒子をゴースト粒子に変更 */
void remove_out_of_bounds(ParticleSystem *ps, double xmin, double xmax,
                          double ymin, double ymax)
{
    double margin = PARTICLE_DISTANCE * (WALL_LAYERS + 1);
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;

        double x = ps->particles[i].pos[0];
        double y = ps->particles[i].pos[1];

        if (x < xmin - margin || x > xmax + margin ||
            y < ymin - margin || y > ymax + margin ||
            isnan(x) || isnan(y)) {
            ps->particles[i].type = GHOST_PARTICLE;
            for (int d = 0; d < DIM; d++) {
                ps->particles[i].vel[d] = 0.0;
                ps->particles[i].acc[d] = 0.0;
            }
        }
    }
}
