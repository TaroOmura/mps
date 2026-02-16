#ifndef CONFIG_H
#define CONFIG_H

/* --- 次元設定 (コンパイル時固定) --- */
#define DIM 3  /* 3次元シミュレーション */

/* --- 粒子種別 (コンパイル時固定) --- */
#define FLUID_PARTICLE    0
#define WALL_PARTICLE     1
#define GHOST_PARTICLE    2

#endif /* CONFIG_H */
