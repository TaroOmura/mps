#include <stdio.h>
#include <math.h>
#include "simulation.h"
#include "operators.h"
#include "pressure_solver.h"
#include "boundary.h"
#include "io.h"
#include "config.h"

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
 *   8. 負圧のクランプ
 *   9. 圧力勾配による速度・位置の修正
 */
void simulation_step(ParticleSystem *ps, NeighborList *nl, int step)
{
    (void)step;
    double re = INFLUENCE_RADIUS;

    /* === 陽的ステップ === */

    /* 加速度の初期化 + 重力の設定 */
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == FLUID_PARTICLE) {
            ps->particles[i].acc[0] = GRAVITY_X;
            ps->particles[i].acc[1] = GRAVITY_Y;
        } else {
            ps->particles[i].acc[0] = 0.0;
            ps->particles[i].acc[1] = 0.0;
        }
    }

    /* 粘性項の計算（accに加算） */
    calc_viscosity_term(ps, nl);

    /* 仮速度の更新 */
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;
        for (int d = 0; d < DIM; d++) {
            ps->particles[i].vel[d] += DT * ps->particles[i].acc[d];
        }
    }

    /* 仮位置の更新 */
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;
        for (int d = 0; d < DIM; d++) {
            ps->particles[i].pos[d] += DT * ps->particles[i].vel[d];
        }
    }

    /* 壁粒子の固定 */
    apply_wall_boundary(ps);

    /* 壁面位置クランプ (陽的ステップ後) */
    clamp_to_walls(ps);

    /* === 陰的ステップ === */

    /* 近傍探索の更新 (仮位置にて) */
    neighbor_search_brute_force(nl, ps, re);

    /* 壁面反発力 */
    apply_wall_repulsion(ps, nl);

    /* 粒子間衝突処理 */
    handle_collision(ps, nl);

    /* 粒子数密度の計算 */
    calc_particle_number_density(ps, nl);

    /* 自由表面判定 */
    judge_free_surface(ps, SURFACE_THRESHOLD);

    /* 圧力ポアソン方程式の求解 */
    solve_pressure(ps, nl);

    /* 負圧のクランプ */
    clamp_negative_pressure(ps);

    /* 圧力勾配による速度・位置の修正 */
    calc_pressure_gradient(ps, nl);

    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != FLUID_PARTICLE) continue;
        for (int d = 0; d < DIM; d++) {
            double du = DT * ps->particles[i].acc[d];
            ps->particles[i].vel[d] += du;
            ps->particles[i].pos[d] += DT * du;
        }
    }

    /* 壁粒子の再固定 */
    apply_wall_boundary(ps);

    /* 壁面位置クランプ (陰的ステップ後) */
    clamp_to_walls(ps);

    /* 領域外粒子の処理 */
    remove_out_of_bounds(ps, DOMAIN_X_MIN, DOMAIN_X_MAX, DOMAIN_Y_MIN, DOMAIN_Y_MAX);
}

/*
 * シミュレーション全体の実行
 */
void simulation_run(ParticleSystem *ps)
{
    double re = INFLUENCE_RADIUS;
    int total_steps = (int)(T_END / DT);

    /* 近傍リストの生成 */
    NeighborList *nl = neighbor_list_create(ps->num, MAX_NEIGHBORS);
    if (!nl) {
        fprintf(stderr, "Error: failed to create neighbor list\n");
        return;
    }

    /* 初期近傍探索 */
    neighbor_search_brute_force(nl, ps, re);

    /* 初期状態の出力 */
    output_csv(ps, 0, "output");
    output_vtk(ps, 0, "output");

    printf("Starting simulation: %d steps, dt = %.2e\n", total_steps, DT);

    for (int step = 1; step <= total_steps; step++) {
        simulation_step(ps, nl, step);

        if (step % OUTPUT_INTERVAL == 0) {
            /* 流体粒子数のカウント */
            int fluid_count = 0;
            for (int i = 0; i < ps->num; i++) {
                if (ps->particles[i].type == FLUID_PARTICLE)
                    fluid_count++;
            }

            printf("Step %6d / %d  (t = %.4f s)  fluid particles: %d\n",
                   step, total_steps, step * DT, fluid_count);
            output_csv(ps, step, "output");
            output_vtk(ps, step, "output");
        }
    }

    printf("Simulation complete.\n");
    neighbor_list_free(nl);
}
