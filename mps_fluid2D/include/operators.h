#ifndef OPERATORS_H
#define OPERATORS_H

#include "particle.h"
#include "neighbor_search.h"

/* 粒子数密度の計算 */
void calc_particle_number_density(ParticleSystem *ps, NeighborList *nl);

/* ラプラシアンモデル: 粘性項の計算 (加速度をaccに加算) */
void calc_viscosity_term(ParticleSystem *ps, NeighborList *nl);

/* 勾配モデル: 圧力勾配による補正加速度をaccに格納 */
void calc_pressure_gradient(ParticleSystem *ps, NeighborList *nl);

/* 自由表面判定 (粒子数密度ベース, 既存手法) */
void judge_free_surface(ParticleSystem *ps, double threshold);

/* 近傍粒子数の計算 (Natsui法用) */
void calc_neighbor_count(ParticleSystem *ps, NeighborList *nl);

/* 自由表面判定 (近傍粒子数ベース, Natsui法) */
void judge_free_surface_by_count(ParticleSystem *ps, double beta);

/* 負圧クランプ: 流体粒子の圧力を max(P, 0) に制限 */
void clamp_pressure(ParticleSystem *ps);

/* 粒子間衝突モデル（越塚 2003）: 陽的ステップ後、近傍探索前に呼ぶ */
void collision(ParticleSystem *ps);

/* ポテンシャル型表面張力 (accに加算) */
void calc_surface_tension(ParticleSystem *ps, NeighborList *nl);

#endif /* OPERATORS_H */
