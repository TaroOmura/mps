#!/usr/bin/env python3
"""
ダムブレイク問題の初期条件ファイルを生成するスクリプト (2D版)

使い方:
    python3 tools/gen_dambreak.py [オプション]

生成されるファイル:
    examples/dambreak/cal.txt        - 計算制御ファイル
    examples/dambreak/params.txt     - パラメータファイル
    examples/dambreak/particles.txt  - 粒子初期条件ファイル
"""
import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="ダムブレイク問題の初期条件生成 (2D)")
    L = 0.146

    # 領域設定
    parser.add_argument("--domain_x", type=float, default=L*4, help="領域x幅 [m]")
    parser.add_argument("--domain_y", type=float, default=L*2, help="領域y高さ [m]")

    # 水柱設定
    parser.add_argument("--water_x", type=float, default=L, help="水柱x幅 [m]")
    parser.add_argument("--water_y", type=float, default=L*2, help="水柱y高さ [m]")

    # 粒子設定
    parser.add_argument("--l0", type=float, default=0.008, help="粒子間距離 [m]")
    parser.add_argument("--wall_layers", type=int, default=4, help="壁粒子の層数")
    parser.add_argument("--dummy_layers", type=int, default=4, help="ダミー粒子の層数")

    # 物性値
    parser.add_argument("--density", type=float, default=1000.0, help="密度 [kg/m^3]")
    parser.add_argument("--viscosity", type=float, default=0, help="動粘性係数 [m^2/s]")
    parser.add_argument("--gravity_y", type=float, default=-9.81, help="重力y成分 [m/s^2]")

    # 時間設定
    parser.add_argument("--dt", type=float, default=5.0e-4, help="時間刻み [s]")
    parser.add_argument("--t_end", type=float, default=2.0, help="終了時刻 [s]")
    parser.add_argument("--output_interval", type=int, default=100, help="出力間隔 [ステップ]")

    # ソルバー設定
    parser.add_argument("--solver_type", type=int, default=1,
                        help="ソルバー種別 (0: CG, 1: ICCG)")
    parser.add_argument("--clamp_negative_pressure", type=int, default=1,
                        help="負圧クランプ (0: 無効, 1: 有効)")

    # 衝突モデル
    parser.add_argument("--restitution_coeff", type=float, default=0.2,
                        help="粒子間衝突の反発係数 e (0: 完全非弾性, 1: 完全弾性)")
    parser.add_argument("--collision_distance_ratio", type=float, default=0.5,
                        help="衝突判定距離の係数 (col_dist = ratio * l0)")

    # PPE定式化
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
    parser.add_argument("--surface_tension_enabled", type=int, default=0,
                        help="表面張力の有効化 (0: 無効, 1: 有効)")
    parser.add_argument("--sigma", type=float, default=0.073,
                        help="表面張力係数 σ [N/m]")
    parser.add_argument("--surface_tension_re_ratio", type=float, default=3.2,
                        help="表面張力影響半径の倍率 (re_st = ratio * l0)")

    # λ計算方法
    parser.add_argument("--use_analytical_lambda", type=int, default=0,
                        help="1: λを解析解で計算, 0: 初期粒子配置から計算 (default: 0)")

    # 出力先
    parser.add_argument("--outdir", type=str, default="examples/dambreak",
                        help="ファイル出力先ディレクトリ")
    parser.add_argument("--output_dir", type=str, default="output/dambreak",
                        help="シミュレーション出力ディレクトリ名")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    l0 = args.l0
    wl = args.wall_layers
    dl = args.dummy_layers
    dx, dy = args.domain_x, args.domain_y
    wx, wy = args.water_x, args.water_y

    # ---- 粒子生成 ----
    particles = []
    n_fluid = 0
    n_wall = 0
    n_dummy = 0

    total_layers = wl + dl
    i_right = int(round(dx / l0))  # 右壁の基準インデックス

    i_min = -(total_layers - 1)
    i_max = i_right + (total_layers - 1)
    j_min = -(total_layers - 1)
    j_max = int(round(dy / l0))

    eps = 1.0e-10

    for i in range(i_min, i_max + 1):
        for j in range(j_min, j_max + 1):
            x = i * l0
            y = j * l0

            ptype = None  # None = 対象外（スキップ）

            # 壁 (境界面上)
            if j <= 0:
                ptype = 1
            if i <= 0:
                ptype = 1
            if i >= i_right:
                ptype = 1

            # ダミー（壁から上塗り）
            if j <= -wl:
                ptype = 3
            if i <= -wl:
                ptype = 3
            if i >= i_right + wl:
                ptype = 3

            # 流体粒子（水柱の内部）
            if ptype is None:
                if (i >= 1 and x <= wx + eps and
                        j >= 1 and y <= wy + eps):
                    ptype = 0

            if ptype is None:
                continue

            particles.append((x, y, 0.0, 0.0, ptype))
            if ptype == 0:
                n_fluid += 1
            elif ptype == 1:
                n_wall += 1
            else:
                n_dummy += 1

    n_total = len(particles)
    print(f"Generated particles: {n_fluid} fluid, {n_wall} wall, "
          f"{n_dummy} dummy, {n_total} total")

    # ---- particles.txt ----
    particle_path = os.path.join(args.outdir, "particles.txt")
    with open(particle_path, "w") as f:
        f.write("# MPS 2D Dam Break - Initial Particle Configuration\n")
        f.write(f"# Domain: {dx} x {dy} m\n")
        f.write(f"# Water column: {wx} x {wy} m\n")
        f.write(f"# Particle distance: {l0} m\n")
        f.write(f"# Fluid: {n_fluid}, Wall: {n_wall}, Dummy: {n_dummy}, Total: {n_total}\n")
        f.write(f"{n_total}\n")
        f.write("# x y vx vy type\n")
        for p in particles:
            f.write(f"{p[0]:.8e} {p[1]:.8e} {p[2]:.8e} {p[3]:.8e} {p[4]}\n")
    print(f"  -> {particle_path}")

    # ---- params.txt ----
    param_path = os.path.join(args.outdir, "params.txt")
    with open(param_path, "w") as f:
        f.write("# MPS 2D Simulation Parameters\n")
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
        f.write(f"use_analytical_lambda {args.use_analytical_lambda}\n")
        f.write(f"cg_max_iter          10000\n")
        f.write(f"cg_tolerance         1.0e-8\n")
        f.write(f"relaxation_coeff     0.2\n")
        f.write(f"clamp_negative_pressure {args.clamp_negative_pressure}\n")
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
        f.write(f"restitution_coeff        {args.restitution_coeff}\n")
        f.write(f"collision_distance_ratio {args.collision_distance_ratio}\n")
        f.write("#\n")
        f.write("# Surface tension\n")
        f.write(f"surface_tension_enabled  {args.surface_tension_enabled}\n")
        f.write(f"surface_tension_coeff    {args.sigma}\n")
        f.write(f"surface_tension_re_ratio {args.surface_tension_re_ratio}\n")
        f.write("#\n")
        f.write("# Domain\n")
        f.write(f"domain_x_min         0.0\n")
        f.write(f"domain_x_max         {dx}\n")
        f.write(f"domain_y_min         0.0\n")
        f.write(f"domain_y_max         {dy}\n")
        f.write("#\n")
        f.write("# Output\n")
        f.write(f"output_dir           {args.output_dir}\n")
    print(f"  -> {param_path}")

    # ---- cal.txt ----
    cal_path = os.path.join(args.outdir, "cal.txt")
    with open(cal_path, "w") as f:
        f.write("# MPS 2D Calculation Control File\n")
        f.write(f"particle_file  particles.txt\n")
        f.write(f"param_file     params.txt\n")
    print(f"  -> {cal_path}")

    print(f"\nRun simulation with:")
    print(f"  ./mps_sim {cal_path}")


if __name__ == "__main__":
    main()
