#include <stdio.h>
#include <math.h>
#include "simulation.h"
#include "operators.h"
#include "pressure_solver.h"
#include "boundary.h"
#include "io.h"
#include "sim_config.h"

/*
 * 1タイムステップの処理 (Semi-implicit法)
 *
 *   [陽的ステップ]
 *   1. 加速度 = 粘性項 + 重力
 *   2. 仮速度 u* = u^n + Δt * acc
 *   3. 仮位置 r* = r^n + Δt * u*
 *
 *   [陰的ステップ]
 *   4. 近傍探索の更新
 *   5. 粒子数密度 n* の計算
 *   6. 自由表面判定
 *   7. 圧力ポアソン方程式の求解
 *   8. 圧力勾配による速度・位置の修正
 */
void simulation_step(ParticleSystem *ps, NeighborList *nl, CellList *cl, int step)
{
    (void)step;
    double re = fmax(g_config->influence_radius_lap, g_config->influence_radius_n);
    double dt = g_config->dt;

    /* === 陽的ステップ === */

    /* 加速度の初期化 + 重力の設定 */
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == FLUID_PARTICLE) {
            for (int d = 0; d < DIM; d++)
                ps->particles[i].acc[d] = g_config->gravity[d];
        } else {
            for (int d = 0; d < DIM; d++)
                ps->particles[i].acc[d] = 0.0;
        }
    }

    /* 粘性項の計算（accに加算） */
    calc_viscosity_term(ps, nl);

    /* 仮速度の更新 */
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;
        for (int d = 0; d < DIM; d++) {
            ps->particles[i].vel[d] += dt * ps->particles[i].acc[d];
        }
    }

    /* 仮位置の更新 */
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;
        for (int d = 0; d < DIM; d++) {
            ps->particles[i].pos[d] += dt * ps->particles[i].vel[d];
        }
    }

    /* 粒子間衝突処理（孤立粒子の重なりによる圧力爆発を防ぐ） */
    collision(ps);

    /* === 陰的ステップ === */

    /* 近傍探索の更新 (仮位置にて) */
    neighbor_search_cell_linked_list(nl, ps, cl, re);

    /* 粒子数密度の計算 */
    calc_particle_number_density(ps, nl);

    /* 自由表面判定 */
    judge_free_surface(ps, g_config->surface_threshold);

    /* 圧力ポアソン方程式の求解 */
    solve_pressure(ps, nl);

    /* 負圧クランプ (clamp_negative_pressure=1 のとき) */
    if (g_config->clamp_negative_pressure)
        clamp_pressure(ps);

    /* 圧力勾配による速度・位置の修正 */
    calc_pressure_gradient(ps, nl);

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;
        for (int d = 0; d < DIM; d++) {
            double du = dt * ps->particles[i].acc[d];
            ps->particles[i].vel[d] += du;
            ps->particles[i].pos[d] += dt * du;
        }
    }

    /* 壁粒子の再固定 */
    apply_wall_boundary(ps);

    /* 領域外粒子の処理 */
    remove_out_of_bounds(ps,
                         g_config->domain_min[0], g_config->domain_max[0],
                         g_config->domain_min[1], g_config->domain_max[1]);
}

/*
 * シミュレーション全体の実行
 */
void simulation_run(ParticleSystem *ps)
{
    double re = fmax(g_config->influence_radius_lap, g_config->influence_radius_n);
    double dt = g_config->dt;
    int total_steps = (int)(g_config->t_end / dt);
    int out_interval = g_config->output_interval;
    const char *out_dir = g_config->output_dir;

    /* 近傍リストの生成 */
    NeighborList *nl = neighbor_list_create(ps->num, g_config->max_neighbors);
    if (!nl) {
        fprintf(stderr, "Error: failed to create neighbor list\n");
        return;
    }

    /* セルリストの生成 (1回だけ確保して毎ステップ再利用) */
    CellList *cl = cell_list_create(ps->num, re,
                                    g_config->domain_min,
                                    g_config->domain_max);
    if (!cl) {
        fprintf(stderr, "Error: failed to create cell list\n");
        neighbor_list_free(nl);
        return;
    }

    /* 初期近傍探索 */
    neighbor_search_cell_linked_list(nl, ps, cl, re);

    /* 初期状態の出力 */
    output_csv(ps, 0, out_dir);
    output_vtk(ps, 0, out_dir);

    printf("Starting simulation (2D): %d steps, dt = %.2e\n", total_steps, dt);

    for (int step = 1; step <= total_steps; step++) {
        simulation_step(ps, nl, cl, step);

        if (step % out_interval == 0) {
            int fluid_count = 0;
            for (int i = 0; i < ps->num; i++) {
                if (ps->particles[i].type == FLUID_PARTICLE)
                    fluid_count++;
            }

            printf("Step %6d / %d  (t = %.4f s)  fluid particles: %d\n",
                   step, total_steps, step * dt, fluid_count);
            output_csv(ps, step, out_dir);
            output_vtk(ps, step, out_dir);
        }
    }

    printf("Simulation complete.\n");
    cell_list_free(cl);
    neighbor_list_free(nl);
}
