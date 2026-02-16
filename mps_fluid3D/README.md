# MPS 3D Fluid Simulation

MPS (Moving Particle Semi-implicit) 法による3次元非圧縮性流体シミュレーションプログラムです。

## 概要

MPS法は粒子法の一種で、流体をラグランジュ的に追跡する粒子で離散化し、ナビエ・ストークス方程式を解きます。本プログラムでは以下の機能を実装しています。

- 3次元非圧縮性流体の数値シミュレーション
- 半陰的圧力解法 (CG法による圧力ポアソン方程式の求解)
- 自由表面の判定 (粒子数密度に基づく)
- 壁境界処理 (壁粒子による反発力モデル)
- 近傍粒子探索
- CSV形式での結果出力

## ディレクトリ構成

```
mps_fluid3D/
├── src/                  # ソースコード (.c)
│   ├── main.c            # メインルーチン・粒子ファイル読み込み
│   ├── simulation.c      # シミュレーションのメインループ
│   ├── particle.c        # 粒子システムの管理
│   ├── neighbor_search.c # 近傍粒子探索
│   ├── kernel.c          # MPS法のカーネル関数
│   ├── operators.c       # 勾配・ラプラシアン等の微分演算子
│   ├── pressure_solver.c # 圧力ポアソン方程式ソルバー (CG法)
│   ├── boundary.c        # 壁境界処理
│   ├── sim_config.c      # 設定の読み込み・管理
│   └── io.c              # ファイル出力
├── include/              # ヘッダファイル (.h)
├── obj/                  # オブジェクトファイル (ビルド時に生成)
├── output/               # シミュレーション結果の出力先
├── tools/
│   └── gen_dambreak.py   # ダムブレイク問題の初期条件生成スクリプト
├── vis/
│   └── visualize.py      # 結果の3D可視化スクリプト
├── cal.txt               # 計算制御ファイル
├── params.txt            # パラメータファイル
├── particles.txt         # 粒子初期条件ファイル
└── Makefile
```

## 必要環境

- GCC (C99以降)
- Make
- Python 3 (初期条件生成・可視化用)
  - NumPy
  - Matplotlib

## ビルド

```bash
make
```

実行ファイル `mps_sim` が生成されます。

```bash
make clean   # ビルド成果物の削除
```

## 使い方

### 1. 初期条件の生成

ダムブレイク問題の初期条件を生成します。

```bash
python3 tools/gen_dambreak.py
```

以下の3ファイルが生成されます。

| ファイル | 内容 |
|---|---|
| `cal.txt` | 計算制御ファイル (粒子ファイル・パラメータファイルのパスを指定) |
| `params.txt` | シミュレーションパラメータ |
| `particles.txt` | 粒子の初期位置・速度・種別 |

主要なオプション:

```bash
python3 tools/gen_dambreak.py \
  --l0 0.025 \          # 粒子間距離 [m]
  --domain_x 0.4 \      # 領域 x幅 [m]
  --domain_y 0.3 \      # 領域 y高さ [m]
  --domain_z 0.2 \      # 領域 z奥行 [m]
  --water_x 0.10 \      # 水柱 x幅 [m]
  --water_y 0.20 \      # 水柱 y高さ [m]
  --water_z 0.20 \      # 水柱 z奥行 [m]
  --dt 5.0e-4 \         # 時間刻み [s]
  --t_end 2.0           # 終了時刻 [s]
```

### 2. シミュレーションの実行

```bash
./mps_sim cal.txt
```

結果は `output/` ディレクトリにCSVファイルとして出力されます。

### 3. 結果の可視化

```bash
python3 vis/visualize.py --input_dir output
```

オプション:

```bash
python3 vis/visualize.py \
  --input_dir output \     # 入力ディレクトリ
  --fluid_only \           # 流体粒子のみ表示
  --save animation.mp4 \   # アニメーション保存 (.mp4 or .gif)
  --interval 100           # フレーム間隔 [ms]
```

## パラメータ

`params.txt` で設定可能な主要パラメータ:

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `particle_distance` | 0.025 | 粒子間距離 [m] |
| `influence_ratio` | 2.1 | 影響半径 / 粒子間距離 |
| `density` | 1000.0 | 密度 [kg/m^3] |
| `viscosity` | 1e-06 | 動粘性係数 [m^2/s] |
| `gravity_y` | -9.81 | 重力加速度 y成分 [m/s^2] |
| `dt` | 0.0005 | 時間刻み [s] |
| `t_end` | 2.0 | 終了時刻 [s] |
| `output_interval` | 100 | 出力間隔 [ステップ] |
| `cg_max_iter` | 10000 | CG法の最大反復回数 |
| `cg_tolerance` | 1.0e-8 | CG法の収束判定閾値 |
| `surface_threshold` | 0.97 | 自由表面の判定閾値 |

## 粒子ファイルの形式

`particles.txt` の形式:

```
# コメント行
粒子数
# x y z vx vy vz type
1.0e-02 2.0e-02 3.0e-02 0.0e+00 0.0e+00 0.0e+00 0
...
```

`type`: 0 = 流体粒子, 1 = 壁粒子

## ライセンス

MIT License
