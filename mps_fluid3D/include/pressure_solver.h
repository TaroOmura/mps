#ifndef PRESSURE_SOLVER_H
#define PRESSURE_SOLVER_H

#include "particle.h"
#include "neighbor_search.h"

/* 圧力ポアソン方程式を構築・求解し、粒子の圧力を更新 */
void solve_pressure(ParticleSystem *ps, NeighborList *nl);

#endif /* PRESSURE_SOLVER_H */
