#!/usr/bin/env python3
"""
無重力下の液滴オシレーション問題の初期条件ファイルを生成するスクリプト (2D版)

正方形ブロックの流体粒子を生成。壁・ダミー粒子は使用しない。
表面張力が有効なとき、正方形→円形への変形・振動が起こる。

使い方:
    python3 tools/gen_droplet.py [オプション]

生成されるファイル:
    examples/droplet/cal.txt        - 計算制御ファイル
    examples/droplet/params.txt     - パラメータファイル
    examples/droplet/particles.txt  - 粒子初期条件ファイル
"""
import argparse
import os
import math


def main():
    parser = argparse.ArgumentParser(description="液滴オシレーション問題の初期条件生成 (2D)")

    # 粒子設定
    parser.add_argument("--l0",  type=float, default=0.002, help="粒子間距離 [m]")
    parser.add_argument("--nx",  type=int,   default=40,    help="x方向の粒子数")
    parser.add_argument("--ny",  type=int,   default=40,    help="y方向の粒子数")

    # 物性値
    parser.add_argument("--density",   type=float, default=1000.0, help="密度 [kg/m^3]")
    parser.add_argument("--viscosity", type=float, default=1.0e-6, help="動粘性係数 [m^2/s]")
    parser.add_argument("--sigma",     type=float, default=0.0728,  help="表面張力係数 σ [N/m]")

    # 時間設定
    parser.add_argument("--dt",              type=float, default=1.0e-4, help="時間刻み [s]")
    parser.add_argument("--t_end",           type=float, default=5.0,   help="終了時刻 [s]")
    parser.add_argument("--output_interval", type=int,   default=100,    help="出力間隔 [ステップ]")

    # ソルバー設定
    parser.add_argument("--solver_type", type=int, default=1,
                        help="ソルバー種別 (0: CG, 1: ICCG)")
    parser.add_argument("--ppe_type", type=int, default=0,
                        help="PPE定式化 (0: 既存密度型, 1: Natsui弱圧縮型)")
    parser.add_argument("--c_ppe", type=float, default=1.01,
                        help="Natsui型PPEの対角係数 c")
    parser.add_argument("--gamma_ppe", type=float, default=0.01,
                        help="Natsui型PPEの密度補正重み γ")

    # 自由表面判定
    parser.add_argument("--surface_detection_method", type=int, default=0,
                        help="自由表面判定法 (0: 粒子数密度, 1: 近傍粒子数)")
    parser.add_argument("--surface_count_threshold", type=float, default=0.85,
                        help="近傍粒子数法の閾値 (method=1 のみ使用)")

    # 表面張力
    parser.add_argument("--surface_tension_re_ratio", type=float, default=3.2,
                        help="表面張力影響半径の倍率 (re_st = ratio * l0)")

    # 出力先
    parser.add_argument("--outdir",     type=str, default="examples/droplet",
                        help="ファイル出力先ディレクトリ")
    parser.add_argument("--output_dir", type=str, default="output/droplet",
                        help="シミュレーション出力ディレクトリ名")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    l0 = args.l0
    nx, ny = args.nx, args.ny

    # 中心が原点になるようにオフセット
    cx = (nx - 1) * l0 / 2.0
    cy = (ny - 1) * l0 / 2.0

    # ---- 粒子生成 ----
    particles = []
    for iy in range(ny):
        for ix in range(nx):
            x = ix * l0 - cx
            y = iy * l0 - cy
            particles.append((x, y, 0.0, 0.0, 0))  # type=0: FLUID

    n_total = len(particles)

    # ---- 参考情報の計算 ----
    side_x = (nx - 1) * l0
    side_y = (ny - 1) * l0
    area = nx * ny * l0 ** 2
    R_eq = math.sqrt(area / math.pi)

    # 毛管安定条件: dt < sqrt(ρ l0³ / (2π σ))
    dt_cap = math.sqrt(args.density * l0 ** 3 / (2.0 * math.pi * args.sigma))

    # Rayleigh 振動周波数 (2D, n=2 楕円モード)
    # ω² = n(n-1)(n+1) σ / (ρ R³)  n=2 → 6σ/(ρ R³)
    omega2 = math.sqrt(6.0 * args.sigma / (args.density * R_eq ** 3))
    T_rayleigh = 2.0 * math.pi / omega2

    # ドメイン: 液滴が動いても出ないよう十分大きく (±5 * R_eq)
    margin = 5.0 * R_eq
    domain_half = max(side_x, side_y) / 2.0 + margin

    print(f"Generated particles: {n_total} fluid (no wall/dummy)")
    print(f"Square block:        {side_x*1e3:.1f} x {side_y*1e3:.1f} mm  "
          f"(nx={nx}, ny={ny}, l0={l0*1e3:.2f} mm)")
    print(f"Equivalent circle R: {R_eq*1e3:.2f} mm")
    print(f"Capillary dt limit:  {dt_cap*1e3:.3f} ms  (dt={args.dt*1e3:.3f} ms)")
    print(f"Rayleigh n=2 mode:   f={omega2/(2*math.pi):.2f} Hz  "
          f"T={T_rayleigh*1e3:.1f} ms  (t_end={args.t_end*1e3:.0f} ms)")

    # ---- particles.txt ----
    particle_path = os.path.join(args.outdir, "particles.txt")
    with open(particle_path, "w") as f:
        f.write("# MPS 2D Droplet Oscillation - Initial Particle Configuration\n")
        f.write(f"# Square block: nx={nx}, ny={ny}, l0={l0} m\n")
        f.write(f"# Equivalent circle radius: {R_eq:.6f} m\n")
        f.write(f"# Rayleigh n=2 period: {T_rayleigh*1e3:.2f} ms\n")
        f.write(f"{n_total}\n")
        f.write("# x y vx vy type\n")
        for p in particles:
            f.write(f"{p[0]:.8e} {p[1]:.8e} {p[2]:.8e} {p[3]:.8e} {p[4]}\n")
    print(f"  -> {particle_path}")

    # ---- params.txt ----
    param_path = os.path.join(args.outdir, "params.txt")
    with open(param_path, "w") as f:
        f.write("# MPS 2D Droplet Oscillation Parameters\n")
        f.write("#\n")
        f.write("# Particle\n")
        f.write(f"particle_distance    {l0}\n")
        f.write(f"influence_ratio_lap  4.0\n")
        f.write(f"influence_ratio_n    2.1\n")
        f.write(f"max_neighbors        512\n")
        f.write(f"wall_layers          0\n")
        f.write(f"dummy_layers         0\n")
        f.write("#\n")
        f.write("# Material\n")
        f.write(f"density              {args.density}\n")
        f.write(f"viscosity            {args.viscosity}\n")
        f.write(f"gravity_x            0.0\n")
        f.write(f"gravity_y            0.0\n")
        f.write("#\n")
        f.write("# Time\n")
        f.write(f"dt                   {args.dt}\n")
        f.write(f"t_end                {args.t_end}\n")
        f.write(f"output_interval      {args.output_interval}\n")
        f.write("#\n")
        f.write("# Pressure solver\n")
        f.write(f"solver_type          {args.solver_type}\n")
        f.write(f"use_analytical_lambda 0\n")
        f.write(f"cg_max_iter          10000\n")
        f.write(f"cg_tolerance         1.0e-8\n")
        f.write(f"relaxation_coeff     0.2\n")
        f.write(f"clamp_negative_pressure 1\n")
        f.write(f"ppe_type             {args.ppe_type}\n")
        f.write(f"c_ppe                {args.c_ppe}\n")
        f.write(f"gamma_ppe            {args.gamma_ppe}\n")
        f.write("#\n")
        f.write("# Free surface\n")
        f.write(f"surface_threshold    0.97\n")
        f.write(f"surface_detection_method {args.surface_detection_method}\n")
        f.write(f"surface_count_threshold  {args.surface_count_threshold}\n")
        f.write("#\n")
        f.write("# Collision model\n")
        f.write(f"restitution_coeff        0.2\n")
        f.write(f"collision_distance_ratio 0.5\n")
        f.write("#\n")
        f.write("# Surface tension\n")
        f.write(f"surface_tension_enabled  1\n")
        f.write(f"surface_tension_coeff    {args.sigma}\n")
        f.write(f"surface_tension_re_ratio {args.surface_tension_re_ratio}\n")
        f.write("#\n")
        f.write("# Domain: large enough for free droplet\n")
        f.write(f"domain_x_min         {-domain_half:.6f}\n")
        f.write(f"domain_x_max         { domain_half:.6f}\n")
        f.write(f"domain_y_min         {-domain_half:.6f}\n")
        f.write(f"domain_y_max         { domain_half:.6f}\n")
        f.write("#\n")
        f.write("# Output\n")
        f.write(f"output_dir           {args.output_dir}\n")
    print(f"  -> {param_path}")

    # ---- cal.txt ----
    cal_path = os.path.join(args.outdir, "cal.txt")
    with open(cal_path, "w") as f:
        f.write("# MPS 2D Droplet Oscillation Control File\n")
        f.write(f"particle_file  particles.txt\n")
        f.write(f"param_file     params.txt\n")
    print(f"  -> {cal_path}")

    print(f"\nRun simulation with:")
    print(f"  ./mps_sim {cal_path}")


if __name__ == "__main__":
    main()
