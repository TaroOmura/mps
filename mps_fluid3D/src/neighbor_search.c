#include <stdlib.h>
#include <string.h>
#include <math.h>
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
 * 影響半径 re 以内の粒子をリストに格納
 */
void neighbor_search_brute_force(NeighborList *nl, ParticleSystem *ps, double re)
{
    int n = ps->num;

    memset(nl->count, 0, n * sizeof(int));

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

            if (r2 < re * re) {
                if (cnt < nl->max_neighbors) {
                    nl->neighbors[i * nl->max_neighbors + cnt] = j;
                    cnt++;
                }
            }
        }
        nl->count[i] = cnt;
    }
}

int neighbor_get(NeighborList *nl, int i, int j)
{
    return nl->neighbors[i * nl->max_neighbors + j];
}
