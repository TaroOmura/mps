#!/usr/bin/env python3
"""
ダムブレイク問題の初期条件ファイルを生成するスクリプト

使い方:
    python3 tools/gen_dambreak.py [オプション]

生成されるファイル:
    cal.txt        - 計算制御ファイル
    params.txt     - パラメータファイル
    particles.txt  - 粒子初期条件ファイル
"""
import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="ダムブレイク問題の初期条件生成")

    # 領域設定
    parser.add_argument("--domain_x", type=float, default=0.4, help="領域x幅 [m]")
    parser.add_argument("--domain_y", type=float, default=0.3, help="領域y高さ [m]")
    parser.add_argument("--domain_z", type=float, default=0.2, help="領域z奥行 [m]")

    # 水柱設定
    parser.add_argument("--water_x", type=float, default=0.10, help="水柱x幅 [m]")
    parser.add_argument("--water_y", type=float, default=0.20, help="水柱y高さ [m]")
    parser.add_argument("--water_z", type=float, default=0.20, help="水柱z奥行 [m]")

    # 粒子設定
    parser.add_argument("--l0", type=float, default=0.025, help="粒子間距離 [m]")
    parser.add_argument("--wall_layers", type=int, default=3, help="壁粒子の層数")

    # 物性値
    parser.add_argument("--density", type=float, default=1000.0, help="密度 [kg/m^3]")
    parser.add_argument("--viscosity", type=float, default=1.0e-6, help="動粘性係数 [m^2/s]")
    parser.add_argument("--gravity_y", type=float, default=-9.81, help="重力y成分 [m/s^2]")

    # 時間設定
    parser.add_argument("--dt", type=float, default=5.0e-4, help="時間刻み [s]")
    parser.add_argument("--t_end", type=float, default=2.0, help="終了時刻 [s]")
    parser.add_argument("--output_interval", type=int, default=100, help="出力間隔 [ステップ]")

    # 出力先
    parser.add_argument("--outdir", type=str, default=".", help="ファイル出力先ディレクトリ")
    parser.add_argument("--output_dir", type=str, default="output", help="シミュレーション出力ディレクトリ名")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    l0 = args.l0
    wl = args.wall_layers
    dx, dy, dz = args.domain_x, args.domain_y, args.domain_z
    wx, wy, wz = args.water_x, args.water_y, args.water_z

    # ---- 粒子生成 ----
    particles = []
    n_fluid = 0
    n_wall = 0

    i_min = -(wl - 1)
    i_max = int(dx / l0) + (wl - 1)
    j_min = -(wl - 1)
    j_max = int(dy / l0)
    k_min = -(wl - 1)
    k_max = int(dz / l0) + (wl - 1)

    eps = 1.0e-10

    for i in range(i_min, i_max + 1):
        for j in range(j_min, j_max + 1):
            for k in range(k_min, k_max + 1):
                x = i * l0
                y = j * l0
                z = k * l0

                is_wall = False

                # 底面
                if j <= 0:
                    is_wall = True
                # 左壁
                if i <= 0 and j > 0:
                    is_wall = True
                # 右壁
                if x >= dx - eps and j > 0:
                    is_wall = True
                # 前壁
                if k <= 0 and j > 0:
                    is_wall = True
                # 後壁
                if z >= dz - eps and j > 0:
                    is_wall = True

                if is_wall:
                    particles.append((x, y, z, 0.0, 0.0, 0.0, 1))
                    n_wall += 1
                else:
                    # 水柱内部
                    if (x > 0.0 and x < wx - eps and
                            y > 0.0 and y < wy - eps and
                            z > 0.0 and z < wz - eps):
                        particles.append((x, y, z, 0.0, 0.0, 0.0, 0))
                        n_fluid += 1

    n_total = len(particles)
    print(f"Generated particles: {n_fluid} fluid, {n_wall} wall, {n_total} total")

    # ---- particles.txt ----
    particle_path = os.path.join(args.outdir, "particles.txt")
    with open(particle_path, "w") as f:
        f.write("# MPS 3D Dam Break - Initial Particle Configuration\n")
        f.write(f"# Domain: {dx} x {dy} x {dz} m\n")
        f.write(f"# Water column: {wx} x {wy} x {wz} m\n")
        f.write(f"# Particle distance: {l0} m\n")
        f.write(f"# Fluid: {n_fluid}, Wall: {n_wall}, Total: {n_total}\n")
        f.write(f"{n_total}\n")
        f.write("# x y z vx vy vz type\n")
        for p in particles:
            f.write(f"{p[0]:.8e} {p[1]:.8e} {p[2]:.8e} "
                    f"{p[3]:.8e} {p[4]:.8e} {p[5]:.8e} {p[6]}\n")
    print(f"  -> {particle_path}")

    # ---- params.txt ----
    wall_repulsion = abs(args.gravity_y) / (2.0 * l0)
    param_path = os.path.join(args.outdir, "params.txt")
    with open(param_path, "w") as f:
        f.write("# MPS 3D Simulation Parameters\n")
        f.write("#\n")
        f.write("# Particle\n")
        f.write(f"particle_distance    {l0}\n")
        f.write(f"influence_ratio      2.1\n")
        f.write(f"max_neighbors        512\n")
        f.write(f"wall_layers          {wl}\n")
        f.write("#\n")
        f.write("# Material\n")
        f.write(f"density              {args.density}\n")
        f.write(f"viscosity            {args.viscosity}\n")
        f.write(f"gravity_x            0.0\n")
        f.write(f"gravity_y            {args.gravity_y}\n")
        f.write(f"gravity_z            0.0\n")
        f.write("#\n")
        f.write("# Time\n")
        f.write(f"dt                   {args.dt}\n")
        f.write(f"t_end                {args.t_end}\n")
        f.write(f"output_interval      {args.output_interval}\n")
        f.write("#\n")
        f.write("# Pressure solver\n")
        f.write(f"cg_max_iter          10000\n")
        f.write(f"cg_tolerance         1.0e-8\n")
        f.write(f"relaxation_coeff     0.2\n")
        f.write("#\n")
        f.write("# Free surface\n")
        f.write(f"surface_threshold    0.97\n")
        f.write("#\n")
        f.write("# Wall\n")
        f.write(f"wall_repulsion_coeff {wall_repulsion:.1f}\n")
        f.write(f"wall_restitution     0.2\n")
        f.write("#\n")
        f.write("# Domain\n")
        f.write(f"domain_x_min         0.0\n")
        f.write(f"domain_x_max         {dx}\n")
        f.write(f"domain_y_min         0.0\n")
        f.write(f"domain_y_max         {dy}\n")
        f.write(f"domain_z_min         0.0\n")
        f.write(f"domain_z_max         {dz}\n")
        f.write("#\n")
        f.write("# Output\n")
        f.write(f"output_dir           {args.output_dir}\n")
    print(f"  -> {param_path}")

    # ---- cal.txt ----
    cal_path = os.path.join(args.outdir, "cal.txt")
    with open(cal_path, "w") as f:
        f.write("# MPS 3D Calculation Control File\n")
        f.write(f"particle_file  particles.txt\n")
        f.write(f"param_file     params.txt\n")
    print(f"  -> {cal_path}")

    print(f"\nRun simulation with:")
    print(f"  ./mps_sim {cal_path}")


if __name__ == "__main__":
    main()
