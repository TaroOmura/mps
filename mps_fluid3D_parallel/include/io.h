#ifndef IO_H
#define IO_H

#include "particle.h"

/* 粒子データをCSV形式で出力 */
/* 列: x, y, z, vx, vy, vz, pressure, type */
void output_csv(const ParticleSystem *ps, int step, const char *output_dir);

/* VTK形式で出力 (ParaView用) */
void output_vtk(const ParticleSystem *ps, int step, const char *output_dir);

#endif /* IO_H */
