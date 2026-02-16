#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "config.h"
#include "sim_config.h"
#include "particle.h"
#include "simulation.h"

/*
 * 粒子初期条件ファイルの読み込み
 *
 * 書式:
 *   # コメント行
 *   粒子数
 *   x y z vx vy vz type
 *   x y z vx vy vz type
 *   ...
 */
static int load_particles(const char *filepath, ParticleSystem *ps)
{
    FILE *fp = fopen(filepath, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open particle file '%s'\n", filepath);
        return -1;
    }

    char line[512];
    int n_particles = 0;
    int header_read = 0;

    /* 粒子数の読み取り (コメント行をスキップ) */
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') continue;

        if (sscanf(p, "%d", &n_particles) == 1) {
            header_read = 1;
            break;
        }
    }

    if (!header_read || n_particles <= 0) {
        fprintf(stderr, "Error: invalid particle count in '%s'\n", filepath);
        fclose(fp);
        return -1;
    }

    int n_fluid = 0, n_wall = 0;
    int loaded = 0;

    while (fgets(line, sizeof(line), fp) && loaded < n_particles) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') continue;

        double x, y, z, vx, vy, vz;
        int type;
        if (sscanf(p, "%lf %lf %lf %lf %lf %lf %d",
                   &x, &y, &z, &vx, &vy, &vz, &type) != 7) {
            fprintf(stderr, "Warning: skipping malformed line: %s", line);
            continue;
        }

        double pos[DIM] = {x, y, z};
        double vel[DIM] = {vx, vy, vz};
        particle_system_add(ps, pos, vel, type);
        loaded++;

        if (type == FLUID_PARTICLE) n_fluid++;
        else if (type == WALL_PARTICLE) n_wall++;
    }

    fclose(fp);

    if (loaded != n_particles) {
        fprintf(stderr, "Warning: expected %d particles but loaded %d\n",
                n_particles, loaded);
    }

    printf("Loaded particles: %d fluid, %d wall, %d total\n",
           n_fluid, n_wall, ps->num);
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <cal_file>\n", argv[0]);
        fprintf(stderr, "  cal_file: calculation file (e.g., cal.txt)\n");
        return 1;
    }

    const char *cal_path = argv[1];

    /* 設定の初期化・読み込み */
    SimConfig config;
    config_set_defaults(&config);

    if (config_load_cal(cal_path, &config) != 0) return 1;
    if (config_load_params(config.param_file, &config) != 0) return 1;

    /* グローバルポインタの設定 */
    g_config = &config;

    printf("=== MPS 3D Fluid Simulation ===\n\n");
    config_print(&config);

    /* 出力ディレクトリの作成 */
    mkdir(config.output_dir, 0755);

    /* 粒子システムの生成 */
    /* 粒子ファイルから行数を事前に取得する代わりに、十分大きな容量を確保 */
    int capacity = 100000;
    ParticleSystem *ps = particle_system_create(capacity);
    if (!ps) {
        fprintf(stderr, "Error: failed to create particle system\n");
        return 1;
    }

    /* 粒子初期条件の読み込み */
    if (load_particles(config.particle_file, ps) != 0) {
        particle_system_free(ps);
        return 1;
    }

    /* 基準粒子数密度 n0 とラプラシアン用パラメータ λ の計算 */
    particle_system_calc_initial_params(ps);

    /* シミュレーション実行 */
    simulation_run(ps);

    /* メモリ解放 */
    particle_system_free(ps);

    return 0;
}
