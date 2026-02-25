#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "sim_config.h"

/* グローバル設定ポインタ */
SimConfig *g_config = NULL;

/* 前後の空白を除去 (ポインタを返す) */
static char *trim(char *s)
{
    while (isspace((unsigned char)*s)) s++;
    if (*s == '\0') return s;
    char *end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) *end-- = '\0';
    return s;
}

void config_set_defaults(SimConfig *config)
{
    memset(config, 0, sizeof(SimConfig));

    config->particle_distance    = 0.025;
    config->influence_ratio_lap  = 2.1;
    config->influence_radius_lap = 2.1 * 0.025;
    config->influence_ratio_n    = 2.1;
    config->influence_radius_n   = 2.1 * 0.025;
    config->max_neighbors        = 256;
    config->wall_layers          = 2;
    config->dummy_layers         = 2;

    config->density   = 1000.0;
    config->viscosity = 1.0e-6;
    config->gravity[0] = 0.0;
    config->gravity[1] = -9.81;

    config->dt              = 5.0e-4;
    config->t_end           = 2.0;
    config->output_interval = 100;

    config->solver_type             = 0;       /* 0: CG, 1: ICCG */
    config->cg_max_iter             = 10000;
    config->cg_tolerance            = 1.0e-8;
    config->relaxation_coeff        = 0.2;
    config->clamp_negative_pressure = 0;

    config->surface_threshold = 0.97;

    config->restitution_coeff        = 0.2;
    config->collision_distance_ratio = 0.5;

    config->domain_min[0] = 0.0;
    config->domain_min[1] = 0.0;
    config->domain_max[0] = 1.0;
    config->domain_max[1] = 0.6;

    config->use_analytical_lambda = 0;

    strncpy(config->output_dir, "output", sizeof(config->output_dir) - 1);
    config->particle_file[0] = '\0';
    config->param_file[0]    = '\0';
}

/*
 * パスの結合: cal.txt のディレクトリを基準にファイルパスを解決
 */
static void resolve_path(const char *cal_path, const char *value,
                         char *out, size_t out_size)
{
    /* 絶対パスならそのまま */
    if (value[0] == '/') {
        strncpy(out, value, out_size - 1);
        out[out_size - 1] = '\0';
        return;
    }

    /* cal_path のディレクトリ部分を取得 */
    char dir[256];
    strncpy(dir, cal_path, sizeof(dir) - 1);
    dir[sizeof(dir) - 1] = '\0';
    char *last_slash = strrchr(dir, '/');
    if (last_slash) {
        *(last_slash + 1) = '\0';
        snprintf(out, out_size, "%s%s", dir, value);
    } else {
        strncpy(out, value, out_size - 1);
        out[out_size - 1] = '\0';
    }
}

/*
 * cal.txt の読み込み
 */
int config_load_cal(const char *cal_path, SimConfig *config)
{
    FILE *fp = fopen(cal_path, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open cal file '%s'\n", cal_path);
        return -1;
    }

    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        char *p = trim(line);
        if (*p == '\0' || *p == '#') continue;

        char key[128], value[256];
        if (sscanf(p, "%127s %255s", key, value) != 2) continue;

        if (strcmp(key, "particle_file") == 0) {
            resolve_path(cal_path, value,
                         config->particle_file, sizeof(config->particle_file));
        } else if (strcmp(key, "param_file") == 0) {
            resolve_path(cal_path, value,
                         config->param_file, sizeof(config->param_file));
        } else {
            fprintf(stderr, "Warning: unknown key in cal file: '%s'\n", key);
        }
    }

    fclose(fp);

    if (config->particle_file[0] == '\0') {
        fprintf(stderr, "Error: particle_file not specified in '%s'\n", cal_path);
        return -1;
    }
    if (config->param_file[0] == '\0') {
        fprintf(stderr, "Error: param_file not specified in '%s'\n", cal_path);
        return -1;
    }

    return 0;
}

/*
 * params.txt の読み込み
 */
int config_load_params(const char *param_path, SimConfig *config)
{
    FILE *fp = fopen(param_path, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open param file '%s'\n", param_path);
        return -1;
    }

    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        char *p = trim(line);
        if (*p == '\0' || *p == '#') continue;

        char key[128], val_str[256];
        if (sscanf(p, "%127s %255s", key, val_str) != 2) continue;

        if (strcmp(key, "particle_distance") == 0) {
            config->particle_distance = atof(val_str);
        } else if (strcmp(key, "influence_ratio_lap") == 0) {
            config->influence_ratio_lap = atof(val_str);
        } else if (strcmp(key, "influence_ratio_n") == 0) {
            config->influence_ratio_n = atof(val_str);
        } else if (strcmp(key, "max_neighbors") == 0) {
            config->max_neighbors = atoi(val_str);
        } else if (strcmp(key, "wall_layers") == 0) {
            config->wall_layers = atoi(val_str);
        } else if (strcmp(key, "dummy_layers") == 0) {
            config->dummy_layers = atoi(val_str);
        } else if (strcmp(key, "density") == 0) {
            config->density = atof(val_str);
        } else if (strcmp(key, "viscosity") == 0) {
            config->viscosity = atof(val_str);
        } else if (strcmp(key, "gravity_x") == 0) {
            config->gravity[0] = atof(val_str);
        } else if (strcmp(key, "gravity_y") == 0) {
            config->gravity[1] = atof(val_str);
        } else if (strcmp(key, "dt") == 0) {
            config->dt = atof(val_str);
        } else if (strcmp(key, "t_end") == 0) {
            config->t_end = atof(val_str);
        } else if (strcmp(key, "output_interval") == 0) {
            config->output_interval = atoi(val_str);
        } else if (strcmp(key, "solver_type") == 0) {
            config->solver_type = atoi(val_str);
        } else if (strcmp(key, "cg_max_iter") == 0) {
            config->cg_max_iter = atoi(val_str);
        } else if (strcmp(key, "cg_tolerance") == 0) {
            config->cg_tolerance = atof(val_str);
        } else if (strcmp(key, "relaxation_coeff") == 0) {
            config->relaxation_coeff = atof(val_str);
        } else if (strcmp(key, "clamp_negative_pressure") == 0) {
            config->clamp_negative_pressure = atoi(val_str);
        } else if (strcmp(key, "surface_threshold") == 0) {
            config->surface_threshold = atof(val_str);
        } else if (strcmp(key, "restitution_coeff") == 0) {
            config->restitution_coeff = atof(val_str);
        } else if (strcmp(key, "collision_distance_ratio") == 0) {
            config->collision_distance_ratio = atof(val_str);
        } else if (strcmp(key, "domain_x_min") == 0) {
            config->domain_min[0] = atof(val_str);
        } else if (strcmp(key, "domain_x_max") == 0) {
            config->domain_max[0] = atof(val_str);
        } else if (strcmp(key, "domain_y_min") == 0) {
            config->domain_min[1] = atof(val_str);
        } else if (strcmp(key, "domain_y_max") == 0) {
            config->domain_max[1] = atof(val_str);
        } else if (strcmp(key, "use_analytical_lambda") == 0) {
            config->use_analytical_lambda = atoi(val_str);
        } else if (strcmp(key, "output_dir") == 0) {
            strncpy(config->output_dir, val_str, sizeof(config->output_dir) - 1);
        } else {
            fprintf(stderr, "Warning: unknown parameter '%s'\n", key);
        }
    }

    fclose(fp);

    /* 影響半径の自動計算 */
    config->influence_radius_lap = config->influence_ratio_lap * config->particle_distance;
    config->influence_radius_n   = config->influence_ratio_n   * config->particle_distance;

    return 0;
}

void config_print(const SimConfig *config)
{
    printf("=== Simulation Configuration ===\n");
    printf("particle_distance:    %.6f m\n", config->particle_distance);
    printf("influence_radius_lap: %.6f m  (ratio = %.2f)  [Laplacian]\n",
           config->influence_radius_lap, config->influence_ratio_lap);
    printf("influence_radius_n:   %.6f m  (ratio = %.2f)  [number density]\n",
           config->influence_radius_n, config->influence_ratio_n);
    printf("max_neighbors:        %d\n", config->max_neighbors);
    printf("wall_layers:          %d\n", config->wall_layers);
    printf("dummy_layers:         %d\n", config->dummy_layers);
    printf("density:              %.1f kg/m^3\n", config->density);
    printf("viscosity:            %.2e m^2/s\n", config->viscosity);
    printf("gravity:              (%.4f, %.4f) m/s^2\n",
           config->gravity[0], config->gravity[1]);
    printf("dt:                   %.2e s\n", config->dt);
    printf("t_end:                %.4f s\n", config->t_end);
    printf("output_interval:      %d steps\n", config->output_interval);
    printf("solver_type:          %s\n", config->solver_type == 1 ? "ICCG" : "CG");
    printf("cg_max_iter:          %d\n", config->cg_max_iter);
    printf("cg_tolerance:         %.2e\n", config->cg_tolerance);
    printf("relaxation_coeff:     %.4f\n", config->relaxation_coeff);
    printf("clamp_negative_pressure: %s\n", config->clamp_negative_pressure ? "ON" : "OFF");
    printf("surface_threshold:    %.4f\n", config->surface_threshold);
    printf("restitution_coeff:         %.4f\n", config->restitution_coeff);
    printf("collision_distance_ratio:  %.4f  (col_dist = %.6f m)\n",
           config->collision_distance_ratio,
           config->collision_distance_ratio * config->particle_distance);
    printf("domain:               [%.3f, %.3f] x [%.3f, %.3f]\n",
           config->domain_min[0], config->domain_max[0],
           config->domain_min[1], config->domain_max[1]);
    printf("use_analytical_lambda: %d  (%s)\n",
           config->use_analytical_lambda,
           config->use_analytical_lambda ? "analytical" : "from initial particles");
    printf("output_dir:           %s\n", config->output_dir);
    printf("particle_file:        %s\n", config->particle_file);
    printf("param_file:           %s\n", config->param_file);
    printf("================================\n\n");
}
