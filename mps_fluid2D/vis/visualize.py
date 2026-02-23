"""
MPS法シミュレーション結果の描画スクリプト
CSVファイルを読み込んで粒子の位置・圧力をアニメーション表示する
"""
import os
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def load_csv(filepath):
    """CSVファイルから粒子データを読み込む"""
    data = np.loadtxt(filepath, delimiter=",", skiprows=1)
    return data


def main():
    parser = argparse.ArgumentParser(description="MPS結果の可視化")
    parser.add_argument("--input_dir", default="output", help="CSVファイルのディレクトリ")
    parser.add_argument("--save", default=None, help="アニメーションを保存するファイルパス (.mp4 or .gif)")
    parser.add_argument("--interval", type=int, default=50, help="フレーム間隔 [ms]")
    parser.add_argument("--dt", type=float, default=5.0e-4, help="シミュレーションの時間刻み [s]")
    parser.add_argument("--fluid_only", action="store_true", help="流体粒子のみを表示する")
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
        ptype = data[:, 5].astype(int)
        all_fluid_pressures.append(data[ptype == 0, 4])
    all_fluid_pressures = np.concatenate(all_fluid_pressures)
    vmin, vmax = np.percentile(all_fluid_pressures, [1, 99])

    # 初期フレームでFigureを構築
    fig, ax = plt.subplots(figsize=(10, 6))
    step0, data0 = frames[0]
    x = data0[:, 0]
    y = data0[:, 1]
    pressure = data0[:, 4]
    ptype = data0[:, 5].astype(int)

    fluid_mask = ptype == 0

    if args.fluid_only:
        wall_scat = None
    else:
        wall_mask = ptype != 0
        wall_scat = ax.scatter(x[wall_mask], y[wall_mask], c="gray", s=2, label="Wall")

    fluid_scat = ax.scatter(x[fluid_mask], y[fluid_mask], c=pressure[fluid_mask],
                            cmap="jet", s=5, vmin=vmin, vmax=vmax, label="Fluid")
    plt.colorbar(fluid_scat, ax=ax, label="Pressure [Pa]")

    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    time0 = step0 * args.dt
    title = ax.set_title(f"MPS Simulation - Step {step0}  (t = {time0:.4f} s)")
    ax.set_aspect("equal")
    ax.legend()
    fig.tight_layout()

    def update(frame_idx):
        step, data = frames[frame_idx]
        x = data[:, 0]
        y = data[:, 1]
        pressure = data[:, 4]
        ptype = data[:, 5].astype(int)

        fluid_mask = ptype == 0

        if wall_scat is not None:
            wall_mask = ptype != 0
            if wall_mask.any():
                wall_scat.set_offsets(np.column_stack([x[wall_mask], y[wall_mask]]))

        if fluid_mask.any():
            fluid_scat.set_offsets(np.column_stack([x[fluid_mask], y[fluid_mask]]))
            fluid_scat.set_array(pressure[fluid_mask])
        else:
            fluid_scat.set_offsets(np.empty((0, 2)))
            fluid_scat.set_array(np.array([]))
        time = step * args.dt
        title.set_text(f"MPS Simulation - Step {step}  (t = {time:.4f} s)")

        if wall_scat is not None:
            return wall_scat, fluid_scat, title
        return fluid_scat, title

    anim = FuncAnimation(fig, update, frames=len(frames),
                         interval=args.interval, blit=True)

    if args.save:
        anim.save(args.save, dpi=150)
        print(f"Saved: {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
