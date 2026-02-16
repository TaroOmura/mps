#ifndef SIMULATION_H
#define SIMULATION_H

#include "particle.h"
#include "neighbor_search.h"

/* シミュレーション全体を実行 */
void simulation_run(ParticleSystem *ps);

/* 1タイムステップの処理 (Semi-implicit法) */
/*   1. 粘性項・外力による仮速度・仮位置の計算 (陽的) */
/*   2. 近傍探索の更新 */
/*   3. 粒子数密度の計算 */
/*   4. 自由表面判定 */
/*   5. 圧力ポアソン方程式の求解 (陰的) */
/*   6. 圧力勾配による速度・位置の修正 */
void simulation_step(ParticleSystem *ps, NeighborList *nl, int step);

#endif /* SIMULATION_H */
