#!/usr/bin/env python3
"""
液滴濡れ問題の初期条件ファイルを生成するスクリプト (2D版)

gen_droplet.py の正方形流体ブロックに、底面の壁粒子・ダミー粒子を追加した構造。
左右・上面は開放。底面のみ壁4層＋ダミー4層。

使い方:
    python3 tools/gen_wetting.py [オプション]

生成されるファイル:
    examples/wetting/cal.txt        - 計算制御ファイル
    examples/wetting/params.txt     - パラメータファイル
    examples/wetting/particles.txt  - 粒子初期条件ファイル
"""
import argparse
import os
import math


def main():
    parser = argparse.ArgumentParser(description="液滴濡れ問題の初期条件生成 (2D)")

    # 粒子設定
    parser.add_argument("--l0",  type=float, default=0.0001, help="粒子間距離 [m]")
    parser.add_argument("--nx",  type=int,   default=21,    help="x方向の粒子数")
    parser.add_argument("--ny",  type=int,   default=21,    help="y方向の粒子数")
    parser.add_argument("--wall_layers",  type=int, default=4, help="壁粒子の層数")
    parser.add_argument("--dummy_layers", type=int, default=4, help="ダミー粒子の層数")
    parser.add_argument("--floor_extra",  type=int, default=20, help="床のX方向の追加マージン [粒子数]")

    # 物性値
    parser.add_argument("--density",   type=float, default=1000.0, help="密度 [kg/m^3]")
    parser.add_argument("--viscosity", type=float, default=1.36e-6,  help="動粘性係数 [m^2/s]")
    parser.add_argument("--gravity_y", type=float, default=-9.81,   help="重力y成分 [m/s^2]")

    # 表面張力
    parser.add_argument("--sigma",     type=float, default=0.0728, help="液液表面張力係数 σ [N/m]")
    parser.add_argument("--surface_tension_re_ratio", type=float, default=3.1,
                        help="表面張力影響半径の倍率 (re_st = ratio * l0)")
    parser.add_argument("--wetting_angle_SL", type=float, default=90.0,
                        help="固液接触角 θ [deg] (0: 完全濡れ, 180: 完全非濡れ)")

    # 時間設定
    parser.add_argument("--dt",              type=float, default=1.0e-5, help="時間刻み [s]")
    parser.add_argument("--t_end",           type=float, default=0.1,   help="終了時刻 [s]")
    parser.add_argument("--output_interval", type=int,   default=100,   help="出力間隔 [ステップ]")

    # ソルバー設定
    parser.add_argument("--solver_type",   type=int, default=1, help="ソルバー種別 (0: CG, 1: ICCG)")
    parser.add_argument("--cmps_gradient", type=int, default=1, help="圧力勾配モデル (0: 標準, 1: CMPS対称型, 2: Oochi対称型)")
    parser.add_argument("--hs_mode",       type=int, default=1, help="HSモード (0: 標準, 1: High order)")
    parser.add_argument("--ppe_type",      type=int, default=0, help="PPE定式化 (0: 密度型, 1: Natsui型)")
    parser.add_argument("--c_ppe",         type=float, default=1.01, help="Natsui型PPEの対角係数 c")
    parser.add_argument("--gamma_ppe",     type=float, default=0.01, help="Natsui型PPEの密度補正重み γ")

    # 自由表面判定
    parser.add_argument("--surface_detection_method", type=int,   default=0,    help="自由表面判定法 (0: 粒子数密度, 1: 近傍粒子数)")
    parser.add_argument("--surface_count_threshold",  type=float, default=0.85, help="近傍粒子数法の閾値")

    # 出力先
    parser.add_argument("--outdir",     type=str, default="examples/wetting", help="ファイル出力先ディレクトリ")
    parser.add_argument("--output_dir", type=str, default="output/wetting",   help="シミュレーション出力ディレクトリ名")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    l0 = args.l0
    nx, ny = args.nx, args.ny
    wl, dl = args.wall_layers, args.dummy_layers

    # x中心オフセット（dropletと同じ）
    cx = (nx - 1) * l0 / 2.0

    # 等価半径の半分を初期ギャップとする（壁から液滴底面までの距離）
    area = nx * ny * l0 ** 2
    R_eq = math.sqrt(area / math.pi)
    gap_l0 = max(1, round(R_eq / (2.0 * l0)))  # l0単位に丸め

    # ---- 粒子生成 ----
    particles = []
    n_fluid = n_wall = n_dummy = 0

    # 流体粒子: 正方形ブロック、底面が y = gap_l0 * l0 になるよう配置
    for iy in range(ny):
        for ix in range(nx):
            x = ix * l0 - cx
            y = (iy + gap_l0) * l0
            particles.append((x, y, 0.0, 0.0, 0))
            n_fluid += 1

    # 壁・ダミー粒子: 底面のみ、x範囲は流体ブロック幅 + 両端 (wl+dl+floor_extra) 層分
    margin = wl + dl + args.floor_extra
    ix_min = -margin
    ix_max = (nx - 1) + margin
    for ix in range(ix_min, ix_max + 1):
        x = ix * l0 - cx
        for jy in range(0, -(wl + dl), -1):  # j=0,-1,...,-(wl+dl-1)
            y = jy * l0
            if jy > -wl:
                ptype = 1  # 壁
                n_wall += 1
            else:
                ptype = 3  # ダミー
                n_dummy += 1
            particles.append((x, y, 0.0, 0.0, ptype))

    n_total = len(particles)

    # ---- 参考情報 ----
    dt_cap = math.sqrt(args.density * l0 ** 3 / (2.0 * math.pi * args.sigma))

    # ドメイン
    wall_x_min = ix_min * l0 - cx
    wall_x_max = ix_max * l0 - cx
    domain_y_min = -(wl + dl) * l0
    domain_y_max = (ny + gap_l0) * l0 + 5.0 * l0  # 流体上端 + 余裕

    print(f"Fluid block: {nx} x {ny}  ({(nx-1)*l0*1e3:.1f} x {(ny-1)*l0*1e3:.1f} mm), "
          f"bottom at y={gap_l0*l0*1e3:.2f} mm  (gap={gap_l0} layers = R_eq/2 ≈ {R_eq/2*1e3:.2f} mm)")
    print(f"Wall: {wl} layers, Dummy: {dl} layers (bottom only)")
    print(f"Equivalent circle R: {R_eq*1e3:.2f} mm")
    print(f"Capillary dt limit:  {dt_cap*1e3:.4f} ms  (dt={args.dt*1e3:.4f} ms)")
    print(f"Wetting angle: {args.wetting_angle_SL} deg")
    print(f"Particles: {n_fluid} fluid, {n_wall} wall, {n_dummy} dummy, {n_total} total")

    # ---- particles.txt ----
    particle_path = os.path.join(args.outdir, "particles.txt")
    with open(particle_path, "w") as f:
        f.write("# MPS 2D Wetting - Initial Particle Configuration\n")
        f.write(f"# Fluid block: nx={nx}, ny={ny}, l0={l0} m\n")
        f.write(f"# Wall layers={wl}, Dummy layers={dl} (bottom only)\n")
        f.write(f"# Wetting angle: {args.wetting_angle_SL} deg\n")
        f.write(f"# Fluid: {n_fluid}, Wall: {n_wall}, Dummy: {n_dummy}, Total: {n_total}\n")
        f.write(f"{n_total}\n")
        f.write("# x y vx vy type\n")
        for p in particles:
            f.write(f"{p[0]:.8e} {p[1]:.8e} {p[2]:.8e} {p[3]:.8e} {p[4]}\n")
    print(f"  -> {particle_path}")

    # ---- params.txt ----
    param_path = os.path.join(args.outdir, "params.txt")
    with open(param_path, "w") as f:
        f.write("# MPS 2D Wetting Parameters\n")
        f.write("#\n")
        f.write("# Particle\n")
        f.write(f"particle_distance    {l0}\n")
        f.write(f"influence_ratio_lap  4.0\n")
        f.write(f"influence_ratio_n    2.1\n")
        f.write(f"max_neighbors        512\n")
        f.write(f"wall_layers          {wl}\n")
        f.write(f"dummy_layers         {dl}\n")
        f.write("#\n")
        f.write("# Material\n")
        f.write(f"density              {args.density}\n")
        f.write(f"viscosity            {args.viscosity}\n")
        f.write(f"gravity_x            0.0\n")
        f.write(f"gravity_y            {args.gravity_y}\n")
        f.write("#\n")
        f.write("# Time\n")
        f.write(f"dt                   {args.dt}\n")
        f.write(f"t_end                {args.t_end}\n")
        f.write(f"output_interval      {args.output_interval}\n")
        f.write("#\n")
        f.write("# Pressure solver\n")
        f.write(f"solver_type          {args.solver_type}\n")
        f.write(f"cmps_gradient        {args.cmps_gradient}\n")
        f.write(f"use_analytical_lambda 0\n")
        f.write(f"hs_mode              {args.hs_mode}\n")
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
        f.write(f"wetting_angle_SL         {args.wetting_angle_SL}\n")
        f.write("#\n")
        f.write("# Domain\n")
        f.write(f"domain_x_min         {wall_x_min:.6f}\n")
        f.write(f"domain_x_max         {wall_x_max:.6f}\n")
        f.write(f"domain_y_min         {domain_y_min:.6f}\n")
        f.write(f"domain_y_max         {domain_y_max:.6f}\n")
        f.write("#\n")
        f.write("# Output\n")
        f.write(f"output_dir           {args.output_dir}\n")
    print(f"  -> {param_path}")

    # ---- cal.txt ----
    cal_path = os.path.join(args.outdir, "cal.txt")
    with open(cal_path, "w") as f:
        f.write("# MPS 2D Wetting Control File\n")
        f.write(f"particle_file  particles.txt\n")
        f.write(f"param_file     params.txt\n")
    print(f"  -> {cal_path}")

    print(f"\nRun simulation with:")
    print(f"  ./mps_sim {cal_path}")


if __name__ == "__main__":
    main()
