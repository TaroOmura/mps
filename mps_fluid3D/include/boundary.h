#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "particle.h"
#include "neighbor_search.h"

/* 壁粒子の境界条件を適用（速度ゼロに固定） */
void apply_wall_boundary(ParticleSystem *ps);


/* 計算領域外に出た粒子の処理 */
void remove_out_of_bounds(ParticleSystem *ps, double xmin, double xmax,
                          double ymin, double ymax,
                          double zmin, double zmax);

#endif /* BOUNDARY_H */
