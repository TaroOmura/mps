"""
MPS法 3Dシミュレーション結果の描画スクリプト
CSVファイルを読み込んで粒子の位置・圧力を3Dアニメーション表示する

※ 3Dデータの詳細な可視化にはParaViewでVTKファイルを開くことを推奨
"""
import os
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib.animation import FuncAnimation


def load_csv(filepath):
    """CSVファイルから粒子データを読み込む"""
    data = np.loadtxt(filepath, delimiter=",", skiprows=1)
    return data


def main():
    parser = argparse.ArgumentParser(description="MPS 3D結果の可視化")
    parser.add_argument("--input_dir", default="../output", help="CSVファイルのディレクトリ")
    parser.add_argument("--save", default=None, help="アニメーションを保存するファイルパス (.mp4 or .gif)")
    parser.add_argument("--interval", type=int, default=100, help="フレーム間隔 [ms]")
    parser.add_argument("--dt", type=float, default=5.0e-4, help="シミュレーションの時間刻み [s]")
    parser.add_argument("--fluid_only", action="store_true", help="流体粒子のみ表示")
    args = parser.parse_args()

    csv_files = sorted(glob.glob(os.path.join(args.input_dir, "*.csv")))
    if not csv_files:
        print(f"CSVファイルが見つかりません: {args.input_dir}")
        return

    # 全フレームのデータを事前読み込み
    frames = []
    for filepath in csv_files:
        filename = os.path.basename(filepath)
        step = int(filename.replace("output_", "").replace(".csv", ""))
        frames.append((step, load_csv(filepath)))
    print(f"{len(frames)} フレームを読み込みました")

    # 圧力の全体レンジを計算（カラーバーの統一用）
    all_fluid_pressures = []
    for step, data in frames:
        ptype = data[:, 7].astype(int)  # type列は8列目 (index 7)
        all_fluid_pressures.append(data[ptype == 0, 6])  # pressure列は7列目 (index 6)
    all_fluid_pressures = np.concatenate(all_fluid_pressures)
    vmin, vmax = np.percentile(all_fluid_pressures, [1, 99])

    # 初期フレームでFigureを構築
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection="3d")

    step0, data0 = frames[0]
    x = data0[:, 0]
    y = data0[:, 1]
    z = data0[:, 2]
    pressure = data0[:, 6]
    ptype = data0[:, 7].astype(int)

    if args.fluid_only:
        fluid_mask = ptype == 0
        scat = ax.scatter(x[fluid_mask], z[fluid_mask], y[fluid_mask],
                          c=pressure[fluid_mask], cmap="jet", s=5,
                          vmin=vmin, vmax=vmax)
    else:
        wall_mask = ptype != 0
        fluid_mask = ptype == 0
        # 壁粒子は半透明の灰色で描画
        ax.scatter(x[wall_mask], z[wall_mask], y[wall_mask],
                   c="gray", s=1, alpha=0.1, label="Wall")
        scat = ax.scatter(x[fluid_mask], z[fluid_mask], y[fluid_mask],
                          c=pressure[fluid_mask], cmap="jet", s=5,
                          vmin=vmin, vmax=vmax, label="Fluid")

    plt.colorbar(scat, ax=ax, label="Pressure [Pa]", shrink=0.6)

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Z [m]")
    ax.set_zlabel("Y [m]")
    time0 = step0 * args.dt
    title = ax.set_title(f"MPS 3D Simulation - Step {step0}  (t = {time0:.4f} s)")
    ax.legend()
    fig.tight_layout()

    def update(frame_idx):
        step, data = frames[frame_idx]
        ax.cla()

        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]
        pressure = data[:, 6]
        ptype = data[:, 7].astype(int)

        fluid_mask = ptype == 0

        if not args.fluid_only:
            wall_mask = ptype != 0
            ax.scatter(x[wall_mask], z[wall_mask], y[wall_mask],
                       c="gray", s=1, alpha=0.1)

        ax.scatter(x[fluid_mask], z[fluid_mask], y[fluid_mask],
                   c=pressure[fluid_mask], cmap="jet", s=5,
                   vmin=vmin, vmax=vmax)

        ax.set_xlabel("X [m]")
        ax.set_ylabel("Z [m]")
        ax.set_zlabel("Y [m]")
        time = step * args.dt
        ax.set_title(f"MPS 3D Simulation - Step {step}  (t = {time:.4f} s)")

    anim = FuncAnimation(fig, update, frames=len(frames),
                         interval=args.interval, blit=False)

    if args.save:
        anim.save(args.save, dpi=150)
        print(f"Saved: {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
