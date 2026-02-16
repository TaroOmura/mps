#include "kernel.h"
#include "config.h"

/* 重み関数の発散を防ぐ最小距離 (粒子間距離の1%) */
#define R_MIN (0.01 * PARTICLE_DISTANCE)

/*
 * MPS標準重み関数
 *   w(r, re) = re/r - 1   (R_MIN < r < re)
 *   w(r, re) = re/R_MIN-1 (r <= R_MIN)  ← 発散防止のクランプ
 *   w(r, re) = 0           (r >= re)
 */
double kernel_weight(double r, double re)
{
    if (r >= re)
        return 0.0;
    if (r < R_MIN)
        r = R_MIN;
    return re / r - 1.0;
}

double kernel_weight_derivative(double r, double re)
{
    if (r >= re)
        return 0.0;
    if (r < R_MIN)
        r = R_MIN;
    return -re / (r * r);
}
