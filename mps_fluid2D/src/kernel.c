#include "kernel.h"
#include "sim_config.h"

/*
 * MPS標準重み関数
 *   w(r, re) = re/r - 1   (r_min < r < re)
 *   w(r, re) = re/r_min-1 (r <= r_min)  ← 発散防止のクランプ
 *   w(r, re) = 0           (r >= re)
 */
double kernel_weight(double r, double re)
{
    double r_min = 0.01 * g_config->particle_distance;
    if (r >= re)
        return 0.0;
    if (r < r_min)
        r = r_min;
    return re / r - 1.0;
}

double kernel_weight_derivative(double r, double re)
{
    double r_min = 0.01 * g_config->particle_distance;
    if (r >= re)
        return 0.0;
    if (r < r_min)
        r = r_min;
    return -re / (r * r);
}
