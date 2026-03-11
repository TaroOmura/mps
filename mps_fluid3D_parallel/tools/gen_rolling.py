#!/usr/bin/env python3
"""
傾いた床での液滴転がりシミュレーション初期条件ファイルを生成するスクリプト (3D版)

傾斜角は重力ベクトルの回転で表現（床は水平、重力をXY平面内で傾ける）。
  gravity_x = g * sin(θ)  … 下流方向 (+x) の加速度成分
  gravity_y = -g * cos(θ) … 床面垂直方向の加速度成分

下流側 (+x) へ床を大きく延長して液滴が転がる様子を観察できるようにする。

使い方:
    python3 tools/gen_rolling.py [オプション]

生成されるファイル:
    examples/rolling/cal.txt        - 計算制御ファイル
    examples/rolling/params.txt     - パラメータファイル
    examples/rolling/particles.txt  - 粒子初期条件ファイル
"""
import argparse
import os
import math


def main():
    parser = argparse.ArgumentParser(description="傾斜面液滴転がりシミュレーション初期条件生成 (3D)")

    # 粒子設定
    parser.add_argument("--l0",  type=float, default=0.0001, help="粒子間距離 [m]")
    parser.add_argument("--nx",  type=int,   default=21,    help="x方向の粒子数")
    parser.add_argument("--ny",  type=int,   default=21,    help="y方向の粒子数")
    parser.add_argument("--nz",  type=int,   default=21,    help="z方向の粒子数")
    parser.add_argument("--wall_layers",  type=int, default=4, help="壁粒子の層数")
    parser.add_argument("--dummy_layers", type=int, default=4, help="ダミー粒子の層数")

    # 傾斜設定
    parser.add_argument("--incline_angle", type=float, default=60.0,
                        help="傾斜角 [deg] (水平面からの角度、重力ベクトルを回転)")

    # 床の非対称延長
    parser.add_argument("--floor_extra_down", type=int, default=60,
                        help="下流側 (+x, 転がり方向) の追加床延長 [粒子数]")
    parser.add_argument("--floor_extra_up",   type=int, default=5,
                        help="上流側 (-x) の追加床延長 [粒子数]")
    parser.add_argument("--floor_extra_side", type=int, default=10,
                        help="側面 (z方向) の追加床延長 [粒子数]")

    # 物性値
    parser.add_argument("--density",           type=float, default=1000.0,  help="密度 [kg/m^3]")
    parser.add_argument("--viscosity",         type=float, default=1.36e-6, help="動粘性係数 [m^2/s]")
    parser.add_argument("--gravity_magnitude", type=float, default=9.81,    help="重力加速度の大きさ [m/s^2]")

    # 表面張力
    parser.add_argument("--sigma",     type=float, default=0.0728, help="液液表面張力係数 σ [N/m]")
    parser.add_argument("--surface_tension_re_ratio", type=float, default=3.1,
                        help="表面張力影響半径の倍率 (re_st = ratio * l0)")
    parser.add_argument("--wetting_angle_SL", type=float, default=120.0,
                        help="固液接触角 θ [deg] (>90: 撥水、液滴が転がりやすい)")

    # 時間設定
    parser.add_argument("--dt",              type=float, default=1.0e-5, help="時間刻み [s]")
    parser.add_argument("--t_end",           type=float, default=0.05,  help="終了時刻 [s]")
    parser.add_argument("--output_interval", type=int,   default=100,   help="出力間隔 [ステップ]")

    # ソルバー設定
    parser.add_argument("--solver_type",   type=int, default=1, help="ソルバー種別 (0: CG, 1: ICCG)")
    parser.add_argument("--cmps_gradient", type=int, default=2, help="圧力勾配モデル (0: 標準, 1: CMPS対称型, 2: Oochi対称型)")
    parser.add_argument("--hs_mode",       type=int, default=1, help="HSモード (0: 標準, 1: High order)")
    parser.add_argument("--ppe_type",      type=int, default=0, help="PPE定式化 (0: 密度型, 1: Natsui型)")
    parser.add_argument("--c_ppe",         type=float, default=1.01, help="Natsui型PPEの対角係数 c")
    parser.add_argument("--gamma_ppe",     type=float, default=0.01, help="Natsui型PPEの密度補正重み γ")

    # 自由表面判定
    parser.add_argument("--surface_detection_method", type=int,   default=0,    help="自由表面判定法 (0: 粒子数密度, 1: 近傍粒子数)")
    parser.add_argument("--surface_count_threshold",  type=float, default=0.85, help="近傍粒子数法の閾値")

    # 出力先
    parser.add_argument("--outdir",     type=str, default="examples/rolling", help="ファイル出力先ディレクトリ")
    parser.add_argument("--output_dir", type=str, default="output/rolling",   help="シミュレーション出力ディレクトリ名")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 傾斜角から重力成分を計算
    angle_rad = math.radians(args.incline_angle)
    gravity_x =  args.gravity_magnitude * math.sin(angle_rad)   # 下流方向 (+x)
    gravity_y = -args.gravity_magnitude * math.cos(angle_rad)   # 床面垂直下向き

    l0 = args.l0
    nx, ny, nz = args.nx, args.ny, args.nz
    wl, dl = args.wall_layers, args.dummy_layers

    # x, z 中心オフセット（液滴は床の上流寄りに配置される）
    cx = (nx - 1) * l0 / 2.0
    cz = (nz - 1) * l0 / 2.0

    # 等価半径の半分を初期ギャップとする
    volume = nx * ny * nz * l0 ** 3
    R_eq = (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0)
    gap_l0 = max(1, round(R_eq / (2.0 * l0)))

    # ---- 粒子生成 ----
    particles = []
    n_fluid = n_wall = n_dummy = 0

    # 流体粒子
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                x = ix * l0 - cx
                y = (iy + gap_l0) * l0
                z = iz * l0 - cz
                particles.append((x, y, z, 0.0, 0.0, 0.0, 0))
                n_fluid += 1

    # 床粒子（非対称延長: 下流方向 (+x) を大きく延長）
    margin_down = wl + dl + args.floor_extra_down
    margin_up   = wl + dl + args.floor_extra_up
    margin_side = wl + dl + args.floor_extra_side
    ix_min = -margin_up
    ix_max = (nx - 1) + margin_down
    iz_min = -margin_side
    iz_max = (nz - 1) + margin_side

    for iz in range(iz_min, iz_max + 1):
        for ix in range(ix_min, ix_max + 1):
            x = ix * l0 - cx
            z = iz * l0 - cz
            for jy in range(0, -(wl + dl), -1):
                y = jy * l0
                if jy > -wl:
                    ptype = 1
                    n_wall += 1
                else:
                    ptype = 3
                    n_dummy += 1
                particles.append((x, y, z, 0.0, 0.0, 0.0, ptype))

    n_total = len(particles)

    # ---- 参考情報 ----
    dt_cap = math.sqrt(args.density * l0 ** 3 / (2.0 * math.pi * args.sigma))
    travel_est = 0.5 * gravity_x * args.t_end ** 2  # 自由加速時の推定移動距離

    wall_x_min = ix_min * l0 - cx
    wall_x_max = ix_max * l0 - cx
    wall_z_min = iz_min * l0 - cz
    wall_z_max = iz_max * l0 - cz
    domain_y_min = -(wl + dl) * l0
    domain_y_max = (ny + gap_l0) * l0 + 10.0 * l0

    print(f"Incline angle: {args.incline_angle} deg  "
          f"(gx={gravity_x:.3f} m/s^2, gy={gravity_y:.3f} m/s^2)")
    print(f"Fluid block: {nx} x {ny} x {nz}  "
          f"({(nx-1)*l0*1e3:.1f} x {(ny-1)*l0*1e3:.1f} x {(nz-1)*l0*1e3:.1f} mm), "
          f"bottom at y={gap_l0*l0*1e3:.2f} mm  (gap={gap_l0} layers)")
    print(f"Floor x: [{wall_x_min*1e3:.2f}, {wall_x_max*1e3:.2f}] mm  "
          f"(up={margin_up} ptcl, down={margin_down} ptcl)")
    print(f"Floor z: [{wall_z_min*1e3:.2f}, {wall_z_max*1e3:.2f}] mm  (side={margin_side} ptcl)")
    print(f"Equivalent sphere R: {R_eq*1e3:.2f} mm")
    print(f"Capillary dt limit:  {dt_cap*1e3:.4f} ms  (dt={args.dt*1e3:.4f} ms)")
    print(f"Wetting angle: {args.wetting_angle_SL} deg")
    print(f"Estimated travel (free accel): {travel_est*1e3:.1f} mm in {args.t_end*1e3:.0f} ms")
    print(f"Particles: {n_fluid} fluid, {n_wall} wall, {n_dummy} dummy, {n_total} total")

    # ---- particles.txt ----
    particle_path = os.path.join(args.outdir, "particles.txt")
    with open(particle_path, "w") as f:
        f.write("# MPS 3D Rolling Droplet - Initial Particle Configuration\n")
        f.write(f"# Fluid block: nx={nx}, ny={ny}, nz={nz}, l0={l0} m\n")
        f.write(f"# Incline angle: {args.incline_angle} deg (gravity rotated in XY plane)\n")
        f.write(f"# Wall layers={wl}, Dummy layers={dl} (bottom only, XZ plane)\n")
        f.write(f"# Floor extra: down={args.floor_extra_down}, up={args.floor_extra_up}, side={args.floor_extra_side}\n")
        f.write(f"# Wetting angle: {args.wetting_angle_SL} deg\n")
        f.write(f"# Fluid: {n_fluid}, Wall: {n_wall}, Dummy: {n_dummy}, Total: {n_total}\n")
        f.write(f"{n_total}\n")
        f.write("# x y z vx vy vz type\n")
        for p in particles:
            f.write(f"{p[0]:.8e} {p[1]:.8e} {p[2]:.8e} "
                    f"{p[3]:.8e} {p[4]:.8e} {p[5]:.8e} {p[6]}\n")
    print(f"  -> {particle_path}")

    # ---- params.txt ----
    param_path = os.path.join(args.outdir, "params.txt")
    with open(param_path, "w") as f:
        f.write("# MPS 3D Rolling Droplet Parameters\n")
        f.write(f"# Incline angle: {args.incline_angle} deg\n")
        f.write("#\n")
        f.write("# Particle\n")
        f.write(f"particle_distance    {l0}\n")
        f.write(f"influence_ratio_lap  3.1\n")
        f.write(f"influence_ratio_n    3.1\n")
        f.write(f"max_neighbors        512\n")
        f.write(f"wall_layers          {wl}\n")
        f.write(f"dummy_layers         {dl}\n")
        f.write("#\n")
        f.write("# Material\n")
        f.write(f"density              {args.density}\n")
        f.write(f"viscosity            {args.viscosity}\n")
        f.write(f"gravity_x            {gravity_x:.8f}\n")
        f.write(f"gravity_y            {gravity_y:.8f}\n")
        f.write(f"gravity_z            0.0\n")
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
        f.write(f"relaxation_coeff     1.0\n")
        f.write(f"clamp_negative_pressure 1\n")
        f.write(f"ppe_type             {args.ppe_type}\n")
        f.write(f"c_ppe                {args.c_ppe}\n")
        f.write(f"gamma_ppe            {args.gamma_ppe}\n")
        f.write("#\n")
        f.write("# Free surface\n")
        f.write(f"surface_threshold    0.85\n")
        f.write(f"surface_detection_method {args.surface_detection_method}\n")
        f.write(f"surface_count_threshold  {args.surface_count_threshold}\n")
        f.write("#\n")
        f.write("# Collision model\n")
        f.write(f"restitution_coeff        0.2\n")
        f.write(f"collision_distance_ratio 0.3\n")
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
        f.write(f"domain_z_min         {wall_z_min:.6f}\n")
        f.write(f"domain_z_max         {wall_z_max:.6f}\n")
        f.write("#\n")
        f.write("# Output\n")
        f.write(f"output_dir           {args.output_dir}\n")
    print(f"  -> {param_path}")

    # ---- cal.txt ----
    cal_path = os.path.join(args.outdir, "cal.txt")
    with open(cal_path, "w") as f:
        f.write("# MPS 3D Rolling Droplet Control File\n")
        f.write(f"particle_file  particles.txt\n")
        f.write(f"param_file     params.txt\n")
    print(f"  -> {cal_path}")

    print(f"\nRun simulation with:")
    print(f"  ./mps_sim {cal_path}")


if __name__ == "__main__":
    main()
