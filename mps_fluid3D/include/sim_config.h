#ifndef SIM_CONFIG_H
#define SIM_CONFIG_H

#include "config.h"

/*
 * 実行時パラメータ構造体
 * cal.txt → params.txt から読み込む
 */
typedef struct {
    /* 粒子パラメータ */
    double particle_distance;   /* 初期粒子間距離 l0 [m] */
    double influence_ratio;     /* 影響半径の倍率 (re = ratio * l0) */
    double influence_radius;    /* 影響半径 re [m] (自動計算) */
    int    max_neighbors;       /* 1粒子あたり最大近傍数 */
    int    wall_layers;         /* 壁粒子の層数 */

    /* 物性値 */
    double density;             /* 流体密度 [kg/m^3] */
    double viscosity;           /* 動粘性係数 [m^2/s] */
    double gravity[DIM];        /* 重力加速度 [m/s^2] */

    /* 時間パラメータ */
    double dt;                  /* 時間刻み [s] */
    double t_end;               /* 終了時刻 [s] */
    int    output_interval;     /* 出力間隔 [ステップ] */

    /* 圧力計算 */
    int    solver_type;         /* 線形ソルバー (0: CG, 1: ICCG) */
    int    cg_max_iter;         /* CG法の最大反復数 */
    double cg_tolerance;        /* CG法の収束判定閾値 */
    double relaxation_coeff;    /* 圧力の緩和係数 */

    /* 自由表面判定 */
    double surface_threshold;   /* n/n0 がこの値未満なら自由表面 */

    /* 壁貫通防止パラメータ */
    double wall_repulsion_coeff; /* 壁面反発力係数 */
    double wall_restitution;     /* 壁面反発時の反発係数 */

    /* 計算領域 */
    double domain_min[DIM];     /* 領域の最小座標 */
    double domain_max[DIM];     /* 領域の最大座標 */

    /* 出力設定 */
    char   output_dir[256];     /* 出力ディレクトリ */

    /* 入力ファイルパス (cal.txt から読み込み) */
    char   particle_file[256];  /* 粒子初期条件ファイル */
    char   param_file[256];     /* パラメータファイル */
} SimConfig;

/* グローバル設定ポインタ (全モジュールから参照可能) */
extern SimConfig *g_config;

/* デフォルト値の設定 */
void config_set_defaults(SimConfig *config);

/* cal.txt の読み込み (particle_file, param_file のパスを取得) */
int config_load_cal(const char *cal_path, SimConfig *config);

/* params.txt の読み込み (計算パラメータを設定) */
int config_load_params(const char *param_path, SimConfig *config);

/* 設定内容の表示 */
void config_print(const SimConfig *config);

#endif /* SIM_CONFIG_H */
