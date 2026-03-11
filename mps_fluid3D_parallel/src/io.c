#include <stdio.h>
#include <string.h>
#include "io.h"
#include "config.h"

/*
 * CSV形式で出力 (3D)
 * 列: x, y, z, vx, vy, vz, pressure, type
 */
void output_csv(const ParticleSystem *ps, int step, const char *output_dir)
{
    char filename[256];
    snprintf(filename, sizeof(filename), "%s/output_%06d.csv", output_dir, step);

    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot open %s\n", filename);
        return;
    }

    fprintf(fp, "x,y,z,vx,vy,vz,pressure,type\n");
    for (int i = 0; i < ps->num; i++) {
        const Particle *p = &ps->particles[i];
        if (p->type == GHOST_PARTICLE) continue;  /* ゴーストは出力しない */
        fprintf(fp, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%d\n",
                p->pos[0], p->pos[1], p->pos[2],
                p->vel[0], p->vel[1], p->vel[2],
                p->pressure, p->type);
    }

    fclose(fp);
}

/*
 * VTK (Legacy) 形式で出力 (ParaView対応, 3D)
 */
void output_vtk(const ParticleSystem *ps, int step, const char *output_dir)
{
    char filename[256];
    snprintf(filename, sizeof(filename), "%s/output_%06d.vtk", output_dir, step);

    /* ゴースト粒子を除いた粒子数をカウント */
    int count = 0;
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type != GHOST_PARTICLE) count++;
    }

    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot open %s\n", filename);
        return;
    }

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "MPS 3D Simulation Step %d\n", step);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    /* 座標 (3D) */
    fprintf(fp, "POINTS %d double\n", count);
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == GHOST_PARTICLE) continue;
        fprintf(fp, "%.8e %.8e %.8e\n",
                ps->particles[i].pos[0], ps->particles[i].pos[1],
                ps->particles[i].pos[2]);
    }

    /* セル（各点が独立セル） */
    fprintf(fp, "CELLS %d %d\n", count, count * 2);
    int idx = 0;
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == GHOST_PARTICLE) continue;
        fprintf(fp, "1 %d\n", idx++);
    }

    fprintf(fp, "CELL_TYPES %d\n", count);
    for (int i = 0; i < count; i++) {
        fprintf(fp, "1\n");  /* VTK_VERTEX */
    }

    /* スカラーデータ */
    fprintf(fp, "POINT_DATA %d\n", count);

    fprintf(fp, "SCALARS pressure double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == GHOST_PARTICLE) continue;
        fprintf(fp, "%.8e\n", ps->particles[i].pressure);
    }

    fprintf(fp, "SCALARS type int 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == GHOST_PARTICLE) continue;
        fprintf(fp, "%d\n", ps->particles[i].type);
    }

    fprintf(fp, "VECTORS velocity double\n");
    for (int i = 0; i < ps->num; i++) {
        if (ps->particles[i].type == GHOST_PARTICLE) continue;
        fprintf(fp, "%.8e %.8e %.8e\n",
                ps->particles[i].vel[0], ps->particles[i].vel[1],
                ps->particles[i].vel[2]);
    }

    fclose(fp);
}
