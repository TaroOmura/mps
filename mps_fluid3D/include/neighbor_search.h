#ifndef NEIGHBOR_SEARCH_H
#define NEIGHBOR_SEARCH_H

#include "particle.h"

/* 近傍リスト構造体 */
typedef struct {
    int *neighbors;   /* 近傍粒子インデックス配列 */
    int *count;       /* 各粒子の近傍粒子数 */
    int  max_neighbors; /* 1粒子あたり最大近傍数 */
} NeighborList;

/* 近傍リストの生成・破棄 */
NeighborList *neighbor_list_create(int num_particles, int max_neighbors);
void          neighbor_list_free(NeighborList *nl);

/* 近傍粒子探索（全探索） */
void neighbor_search_brute_force(NeighborList *nl, ParticleSystem *ps, double re);

/* i番目の粒子のj番目の近傍インデックスを取得 */
int neighbor_get(NeighborList *nl, int i, int j);

#endif /* NEIGHBOR_SEARCH_H */
