#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "particle.h"
#include "simulation.h"

/*
 * ダムブレイク問題の初期粒子配置
 *
 *     0.0     WATER_X              DOMAIN_X_MAX
 * 0.6  +---+---------------------------+
 *      | W |                           |
 * 0.5  | A |   (air)                   |
 *      | T |                           |
 *      | E |                           |
 *      | R |                           |
 *      +---+                           |
 *      |                               |
 * 0.0  +-------------------------------+
 *
 * 壁粒子: 底面・左壁・右壁に WALL_LAYERS 層
 * 流体粒子: 水柱内部
 */
static void setup_dam_break(ParticleSystem *ps)
{
    double l0 = PARTICLE_DISTANCE;
    double vel[DIM] = {0.0, 0.0};
    int n_fluid = 0, n_wall = 0;

    /*
     * 全粒子を格子点 (i*l0, j*l0) に配置
     * i, j の範囲は壁の外側層も含む
     */
    int i_min = -(WALL_LAYERS - 1);
    int i_max = (int)(DOMAIN_X_MAX / l0) + (WALL_LAYERS - 1);
    int j_min = -(WALL_LAYERS - 1);
    int j_max = (int)(DOMAIN_Y_MAX / l0);

    for (int i = i_min; i <= i_max; i++) {
        for (int j = j_min; j <= j_max; j++) {
            double x = i * l0;
            double y = j * l0;
            double pos[DIM] = {x, y};

            /* 壁領域の判定: 底面(y<=0), 左壁(x<=0), 右壁(x>=DOMAIN_X_MAX) */
            int is_wall = 0;

            /* 底面の壁 */
            if (j <= 0) {
                is_wall = 1;
            }
            /* 左壁 */
            if (i <= 0 && j > 0) {
                is_wall = 1;
            }
            /* 右壁 */
            if (x >= DOMAIN_X_MAX - 1.0e-10 && j > 0) {
                is_wall = 1;
            }

            if (is_wall) {
                particle_system_add(ps, pos, vel, WALL_PARTICLE);
                n_wall++;
            } else {
                /* 流体領域: 水柱内部 */
                if (x > 0.0 && x < WATER_X - 1.0e-10 &&
                    y > 0.0 && y < WATER_Y - 1.0e-10) {
                    particle_system_add(ps, pos, vel, FLUID_PARTICLE);
                    n_fluid++;
                }
                /* それ以外（空気領域）は粒子を配置しない */
            }
        }
    }

    printf("Dam break setup: %d fluid, %d wall, %d total\n",
           n_fluid, n_wall, ps->num);
}

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    printf("=== MPS Fluid Simulation ===\n");
    printf("Particle distance: %.4f m\n", PARTICLE_DISTANCE);
    printf("Influence radius:  %.4f m\n", INFLUENCE_RADIUS);
    printf("Time step:         %.2e s\n", DT);
    printf("End time:          %.2f s\n", T_END);
    printf("\n");

    /* 粒子システムの生成 */
    ParticleSystem *ps = particle_system_create(MAX_PARTICLES);
    if (!ps) {
        fprintf(stderr, "Error: failed to create particle system\n");
        return 1;
    }

    /* ダムブレイク問題の初期配置 */
    setup_dam_break(ps);

    /* 基準粒子数密度 n0 とラプラシアン用パラメータ λ の計算 */
    particle_system_calc_initial_params(ps);

    /* シミュレーション実行 */
    simulation_run(ps);

    /* メモリ解放 */
    particle_system_free(ps);

    return 0;
}
