#ifndef NEIGHBOR_SEARCH_H
#define NEIGHBOR_SEARCH_H

#include "particle.h"
#include "config.h"

/* 近傍リスト構造体 */
typedef struct {
    int *neighbors;     /* 近傍粒子インデックス配列 [num_particles × max_neighbors] */
    int *count;         /* 各粒子の近傍粒子数 */
    int  max_neighbors; /* 1粒子あたり最大近傍数 */
} NeighborList;

/*
 * セルリンクリスト構造体
 *
 * 計算領域をセルサイズ cell_size(= re)で分割したグリッドを管理する。
 * head[ci] : セル ci の先頭粒子インデックス (-1 = 空)
 * next[i]  : 粒子 i と同一セル内の次粒子インデックス (-1 = 末尾)
 */
typedef struct {
    int   *head;         /* head[cell_idx]: セル内先頭粒子, -1 = 空 */
    int   *next;         /* next[i]: 同セル内の次粒子, -1 = 末尾 */
    int    nx;           /* x方向のセル数 */
    int    ny;           /* y方向のセル数 */
    int    nz;           /* z方向のセル数 */
    int    total_cells;  /* nx * ny * nz */
    double cell_size;    /* セルサイズ [m] (= re) */
    double origin[DIM];  /* グリッド原点座標 */
} CellList;

/* 近傍リストの生成・破棄 */
NeighborList *neighbor_list_create(int num_particles, int max_neighbors);
void          neighbor_list_free(NeighborList *nl);

/* セルリストの生成・破棄 */
CellList *cell_list_create(int num_particles, double re,
                           const double domain_min[DIM],
                           const double domain_max[DIM]);
void      cell_list_free(CellList *cl);

/* 近傍粒子探索（全探索 O(N²)） */
void neighbor_search_brute_force(NeighborList *nl, ParticleSystem *ps, double re);

/* 近傍粒子探索（セルリンクリスト O(N)） */
void neighbor_search_cell_linked_list(NeighborList *nl, ParticleSystem *ps,
                                      CellList *cl, double re);

/* i番目の粒子のj番目の近傍インデックスを取得 */
int neighbor_get(NeighborList *nl, int i, int j);

#endif /* NEIGHBOR_SEARCH_H */
