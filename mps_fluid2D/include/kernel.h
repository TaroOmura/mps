#ifndef KERNEL_H
#define KERNEL_H

/* MPS重み関数 w(r, re) */
double kernel_weight(double r, double re);

/* 重み関数の勾配に使う係数計算 */
double kernel_weight_derivative(double r, double re);

#endif /* KERNEL_H */
