#include "kernel.h"
#include "sim_config.h"

/*
 * MPS標準重み関数
 *   w(r, re) = re/r - 1   (R_MIN < r < re)
 *   w(r, re) = re/R_MIN-1 (r <= R_MIN)  ← 発散防止のクランプ
 *   w(r, re) = 0           (r >= re)
 *
 * R_MIN = 0.01 * particle_distance
 */
double kernel_weight(double r, double re)
{
    if (r >= re)
        return 0.0;
    double r_min = 0.01 * g_config->particle_distance;
    if (r < r_min)
        r = r_min;
    return re / r - 1.0;
}

double kernel_weight_derivative(double r, double re)
{
    if (r >= re)
        return 0.0;
    double r_min = 0.01 * g_config->particle_distance;
    if (r < r_min)
        r = r_min;
    return -re / (r * r);
}
