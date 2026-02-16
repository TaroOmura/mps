#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "particle.h"
#include "neighbor_search.h"

/* 壁粒子の境界条件を適用（速度ゼロに固定） */
void apply_wall_boundary(ParticleSystem *ps);

/* 自由表面の圧力をゼロに設定 */
void apply_surface_pressure(ParticleSystem *ps);

/*
 * 粒子間衝突処理 (Koshizuka & Oka, 1996)
 * 過度に接近した粒子対の速度を補正して重なりを防止
 */
void handle_collision(ParticleSystem *ps, NeighborList *nl);

/* 壁面反発力: 壁粒子から近距離の流体粒子に反発加速度を適用 */
void apply_wall_repulsion(ParticleSystem *ps, NeighborList *nl);

/* 壁面位置クランプ: 壁ジオメトリを貫通した粒子を押し戻す */
void clamp_to_walls(ParticleSystem *ps);

/* 計算領域外に出た粒子の処理 */
void remove_out_of_bounds(ParticleSystem *ps, double xmin, double xmax,
                          double ymin, double ymax);

#endif /* BOUNDARY_H */
