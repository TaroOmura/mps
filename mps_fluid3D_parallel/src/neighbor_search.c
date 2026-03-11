/* OpenMP parallelization added */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "neighbor_search.h"
#include "config.h"

NeighborList *neighbor_list_create(int num_particles, int max_neighbors)
{
    NeighborList *nl = (NeighborList *)malloc(sizeof(NeighborList));
    if (!nl) return NULL;

    nl->neighbors = (int *)malloc(num_particles * max_neighbors * sizeof(int));
    nl->count = (int *)calloc(num_particles, sizeof(int));
    nl->max_neighbors = max_neighbors;

    if (!nl->neighbors || !nl->count) {
        free(nl->neighbors);
        free(nl->count);
        free(nl);
        return NULL;
    }
    return nl;
}

void neighbor_list_free(NeighborList *nl)
{
    if (!nl) return;
    free(nl->neighbors);
    free(nl->count);
    free(nl);
}

/*
 * 全探索による近傍粒子探索 O(N^2)
 *
 * OpenMP並列化:
 *   外側ループ（粒子i）を並列化。各スレッドは独立した粒子iを担当し、
 *   nl->count[i] および nl->neighbors[i * max_neighbors + ...] への書き込みは
 *   スレッド間で衝突しない（各iが独立したメモリ領域を持つ）。
 */
void neighbor_search_brute_force(NeighborList *nl, ParticleSystem *ps, double re)
{
    int n = ps->num;
    int max_nb = nl->max_neighbors;
    double re2 = re * re;

    memset(nl->count, 0, n * sizeof(int));

    #pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type == GHOST_PARTICLE) continue;

        int cnt = 0;
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            if (ps->particles[j].type == GHOST_PARTICLE) continue;

            double r2 = 0.0;
            for (int d = 0; d < DIM; d++) {
                double diff = ps->particles[j].pos[d] - ps->particles[i].pos[d];
                r2 += diff * diff;
            }

            if (r2 < re2) {
                if (cnt >= max_nb) {
                    fprintf(stderr,
                            "Error: neighbor count exceeded max_neighbors (%d) for particle %d. "
                            "Increase max_neighbors in params.txt.\n",
                            max_nb, i);
                    exit(EXIT_FAILURE);
                }
                nl->neighbors[i * max_nb + cnt] = j;
                cnt++;
            }
        }
        nl->count[i] = cnt;
    }
}

int neighbor_get(NeighborList *nl, int i, int j)
{
    return nl->neighbors[i * nl->max_neighbors + j];
}

/* ------------------------------------------------------------------ */
/* セルリンクリスト                                                      */
/* ------------------------------------------------------------------ */

CellList *cell_list_create(int num_particles, double re,
                           const double domain_min[DIM],
                           const double domain_max[DIM])
{
    CellList *cl = (CellList *)malloc(sizeof(CellList));
    if (!cl) return NULL;

    cl->cell_size = re;

    double margin = 4.0 * re;
    for (int d = 0; d < DIM; d++)
        cl->origin[d] = domain_min[d] - margin;

    double span_x = (domain_max[0] + margin) - cl->origin[0];
    double span_y = (domain_max[1] + margin) - cl->origin[1];
    double span_z = (domain_max[2] + margin) - cl->origin[2];

    cl->nx = (int)(span_x / re) + 2;
    cl->ny = (int)(span_y / re) + 2;
    cl->nz = (int)(span_z / re) + 2;
    cl->total_cells = cl->nx * cl->ny * cl->nz;

    cl->head = (int *)malloc(cl->total_cells * sizeof(int));
    cl->next = (int *)malloc(num_particles  * sizeof(int));

    if (!cl->head || !cl->next) {
        cell_list_free(cl);
        return NULL;
    }
    return cl;
}

void cell_list_free(CellList *cl)
{
    if (!cl) return;
    free(cl->head);
    free(cl->next);
    free(cl);
}

/*
 * neighbor_search_cell_linked_list (OpenMP並列版)
 *
 * アルゴリズム:
 *   フェーズ1 (セルリスト構築): cl->head[] および cl->next[] の構築は
 *     複数スレッドが同一セルの head に書き込む可能性があるため、
 *     フェーズ1は逐次処理のまま維持する。
 *     （原子演算でも構築できるが、フェーズ1のコストはフェーズ2に比べて小さいため
 *       逐次処理で十分）
 *
 *   フェーズ2 (近傍探索): 各粒子iは独立したメモリ領域に書き込むため並列化可能。
 *     nl->count[i] と nl->neighbors[i * max_neighbors + ...] はスレッド間で衝突しない。
 *     cl->head[] と cl->next[] は読み取り専用なので安全。
 */
void neighbor_search_cell_linked_list(NeighborList *nl, ParticleSystem *ps,
                                      CellList *cl, double re)
{
    int    n      = ps->num;
    double re2    = re * re;
    int    nx     = cl->nx;
    int    ny     = cl->ny;
    int    nz     = cl->nz;
    int    max_nb = nl->max_neighbors;
    double cs_inv = 1.0 / cl->cell_size;

    /* 近傍数リセット */
    memset(nl->count, 0, n * sizeof(int));

    /* ---- フェーズ1: セルリストの構築 (逐次) ----
     * head[] への prepend は競合が発生するため逐次で実行する。
     */
    memset(cl->head, -1, cl->total_cells * sizeof(int));

    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type == GHOST_PARTICLE) {
            cl->next[i] = -1;
            continue;
        }

        int ix = (int)((ps->particles[i].pos[0] - cl->origin[0]) * cs_inv);
        int iy = (int)((ps->particles[i].pos[1] - cl->origin[1]) * cs_inv);
        int iz = (int)((ps->particles[i].pos[2] - cl->origin[2]) * cs_inv);

        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz) {
            cl->next[i] = -1;
            continue;
        }

        int ci       = (iz * ny + iy) * nx + ix;
        cl->next[i]  = cl->head[ci];
        cl->head[ci] = i;
    }

    /* ---- フェーズ2: 近傍探索 (並列) ----
     * 各粒子 i の書き込み先 (nl->count[i], nl->neighbors[i*max_nb..]) は独立。
     * cl->head[] と cl->next[] は読み取り専用 → データ競合なし。
     */
    #pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < n; i++) {
        if (ps->particles[i].type == GHOST_PARTICLE) continue;

        int ix = (int)((ps->particles[i].pos[0] - cl->origin[0]) * cs_inv);
        int iy = (int)((ps->particles[i].pos[1] - cl->origin[1]) * cs_inv);
        int iz = (int)((ps->particles[i].pos[2] - cl->origin[2]) * cs_inv);

        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz) continue;

        int cnt = 0;

        int cx_min = (ix > 0)      ? ix - 1 : 0;
        int cx_max = (ix < nx - 1) ? ix + 1 : nx - 1;
        int cy_min = (iy > 0)      ? iy - 1 : 0;
        int cy_max = (iy < ny - 1) ? iy + 1 : ny - 1;
        int cz_min = (iz > 0)      ? iz - 1 : 0;
        int cz_max = (iz < nz - 1) ? iz + 1 : nz - 1;

        for (int cz = cz_min; cz <= cz_max; cz++) {
            for (int cy = cy_min; cy <= cy_max; cy++) {
                for (int cx = cx_min; cx <= cx_max; cx++) {
                    int j = cl->head[(cz * ny + cy) * nx + cx];
                    while (j >= 0) {
                        if (j != i && ps->particles[j].type != GHOST_PARTICLE) {
                            double r2 = 0.0;
                            for (int d = 0; d < DIM; d++) {
                                double diff = ps->particles[j].pos[d]
                                            - ps->particles[i].pos[d];
                                r2 += diff * diff;
                            }
                            if (r2 < re2) {
                                if (cnt >= max_nb) {
                                    fprintf(stderr,
                                        "Error: neighbor count exceeded"
                                        " max_neighbors (%d) for particle %d."
                                        " Increase max_neighbors in params.txt.\n",
                                        max_nb, i);
                                    exit(EXIT_FAILURE);
                                }
                                nl->neighbors[i * max_nb + cnt] = j;
                                cnt++;
                            }
                        }
                        j = cl->next[j];
                    }
                }
            }
        }
        nl->count[i] = cnt;
    }
}
