#ifndef PARTICLE_H
#define PARTICLE_H

#include "config.h"

/* 粒子データ構造体 */
typedef struct {
    double pos[DIM];    /* 位置 */
    double vel[DIM];    /* 速度 */
    double acc[DIM];    /* 加速度 */
    double pressure;    /* 圧力 */
    double n;           /* 粒子数密度 */
    int    neighbor_count; /* 近傍粒子数 Ni (Natsui法の自由表面判定用) */
    int    type;        /* 粒子種別 (FLUID / WALL / GHOST) */
    int    on_surface;  /* 自由表面フラグ */
} Particle;

/* 粒子群管理構造体 */
typedef struct {
    Particle *particles;   /* 粒子配列 */
    int       num;         /* 粒子総数 */
    int       capacity;    /* 配列の確保容量 */
    double    n0;          /* 基準粒子数密度 */
    double    lambda;      /* ラプラシアンモデル用パラメータ */
    double    C_LL;        /* 表面張力係数 (初期化時に計算・保存) */
    int       n0_count;    /* 基準近傍粒子数 N0 (Natsui法の自由表面判定用) */
} ParticleSystem;

/* 粒子システムの生成・破棄 */
ParticleSystem *particle_system_create(int capacity);
void            particle_system_free(ParticleSystem *ps);

/* 粒子の追加 */
int particle_system_add(ParticleSystem *ps, double pos[], double vel[], int type);

/* 基準値(n0, lambda)の計算 */
void particle_system_calc_initial_params(ParticleSystem *ps);

#endif /* PARTICLE_H */
