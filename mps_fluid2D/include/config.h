#ifndef CONFIG_H
#define CONFIG_H

/* --- 次元設定 --- */
#define DIM 2  /* 2次元シミュレーション */

/* --- 粒子パラメータ --- */
#define PARTICLE_DISTANCE  0.025   /* 初期粒子間距離 l0 [m] */
#define INFLUENCE_RADIUS   (2.1 * PARTICLE_DISTANCE)  /* 影響半径 re */
#define MAX_PARTICLES      10000
#define MAX_NEIGHBORS      256     /* 1粒子あたり最大近傍数 */
#define WALL_LAYERS        3       /* 壁粒子の層数 */

/* --- 物性値 --- */
#define DENSITY    1000.0  /* 流体密度 [kg/m^3] */
#define VISCOSITY  1.0e-6  /* 動粘性係数 [m^2/s] */
#define GRAVITY_X  0.0     /* 重力加速度 x成分 [m/s^2] */
#define GRAVITY_Y -9.81    /* 重力加速度 y成分 [m/s^2] */

/* --- 時間パラメータ --- */
#define DT              5.0e-4  /* 時間刻み [s] */
#define T_END           2.0     /* 終了時刻 [s] */
#define OUTPUT_INTERVAL 100     /* 出力間隔 [ステップ] */

/* --- 圧力計算 --- */
#define CG_MAX_ITER    10000   /* CG法の最大反復数 */
#define CG_TOLERANCE   1.0e-8  /* CG法の収束判定閾値 */
#define RELAXATION_COEFF 0.2   /* 圧力の緩和係数 */

/* --- 自由表面判定 --- */
#define SURFACE_THRESHOLD 0.97  /* n/n0 がこの値未満なら自由表面 */

/* --- 粒子種別 --- */
#define FLUID_PARTICLE    0
#define WALL_PARTICLE     1
#define GHOST_PARTICLE    2

/* --- 壁貫通防止パラメータ --- */
#define WALL_REPULSION_COEFF 196.2   /* 壁面反発力係数 ≒ 0.5*|g|/l0 */
#define WALL_RESTITUTION     0.2     /* 壁面反発時の反発係数 */

/* --- ダムブレイク問題の領域設定 --- */
#define DOMAIN_X_MIN  0.0
#define DOMAIN_X_MAX  1.0
#define DOMAIN_Y_MIN  0.0
#define DOMAIN_Y_MAX  0.6
#define WATER_X       0.25   /* 水柱の幅 [m] */
#define WATER_Y       0.50   /* 水柱の高さ [m] */

#endif /* CONFIG_H */
