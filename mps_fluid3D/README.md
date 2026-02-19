# MPS 3D Fluid Simulation

MPS (Moving Particle Semi-implicit) 法による3次元非圧縮性流体シミュレーションプログラムです。

## 概要

MPS法は粒子法の一種で、流体をラグランジュ的に追跡する粒子で離散化し、ナビエ・ストークス方程式を解きます。本プログラムでは以下の機能を実装しています。

- 3次元非圧縮性流体の数値シミュレーション
- 半陰的圧力解法 (CG法 / ICCG法を選択可能)
- 自由表面の判定 (粒子数密度に基づく)
- 壁境界処理 (壁粒子による反発力モデル)
- 近傍粒子探索
- CSV形式 / VTK形式での結果出力

## ディレクトリ構成

```
mps_fluid3D/
├── src/                         # ソースコード (.c)
│   ├── main.c                   # メインルーチン・粒子ファイル読み込み
│   ├── simulation.c             # シミュレーションのメインループ
│   ├── particle.c               # 粒子システムの管理
│   ├── neighbor_search.c        # 近傍粒子探索
│   ├── kernel.c                 # MPS法のカーネル関数
│   ├── operators.c              # 勾配・ラプラシアン等の微分演算子
│   ├── pressure_solver.c        # 圧力ポアソン方程式ソルバー (CG法 / ICCG法)
│   ├── boundary.c               # 壁境界処理
│   ├── sim_config.c             # 設定の読み込み・管理
│   └── io.c                     # ファイル出力 (CSV / VTK)
├── include/                     # ヘッダファイル (.h)
├── obj/                         # オブジェクトファイル (ビルド時に生成, git管理外)
├── output/                      # シミュレーション結果の出力先 (git管理外)
├── tools/
│   └── gen_dambreak.py          # ダムブレイク問題の初期条件生成スクリプト
├── vis/
│   └── visualize.py             # 結果の3D可視化スクリプト
├── examples/
│   └── dambreak/                # ダムブレイク問題の入力ファイル
│       ├── cal.txt              # 計算制御ファイル
│       ├── params.txt           # パラメータファイル
│       └── particles.txt        # 粒子初期条件ファイル
├── Makefile
├── .gitignore
└── README.md
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

デフォルトでは `examples/dambreak/` に以下の3ファイルが生成されます。

| ファイル | 内容 |
|---|---|
| `examples/dambreak/cal.txt` | 計算制御ファイル (粒子ファイル・パラメータファイルのパスを指定) |
| `examples/dambreak/params.txt` | シミュレーションパラメータ |
| `examples/dambreak/particles.txt` | 粒子の初期位置・速度・種別 |

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
./mps_sim examples/dambreak/cal.txt
```

または

```bash
make run
```

結果は `output/` ディレクトリにCSVファイルとVTKファイルとして出力されます。
VTKファイルは [ParaView](https://www.paraview.org/) で開いて3Dデータを詳細に可視化できます。

### 3. 結果の可視化

```bash
python3 vis/visualize.py
```

オプション:

```bash
python3 vis/visualize.py \
  --input_dir output \     # 入力ディレクトリ (デフォルト: output)
  --fluid_only \           # 流体粒子のみ表示
  --save animation.mp4 \   # アニメーション保存 (.mp4 or .gif)
  --interval 100           # フレーム間隔 [ms]
```

## パラメータ

`params.txt` で設定可能なパラメータ:

### 粒子設定

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `particle_distance` | 0.025 | 粒子間距離 [m] |
| `influence_ratio` | 2.1 | 影響半径 / 粒子間距離 |
| `max_neighbors` | 512 | 1粒子あたりの最大近傍数 |
| `wall_layers` | 3 | 壁粒子の層数 |

### 物性値

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `density` | 1000.0 | 密度 [kg/m^3] |
| `viscosity` | 1e-06 | 動粘性係数 [m^2/s] |
| `gravity_x` | 0.0 | 重力加速度 x成分 [m/s^2] |
| `gravity_y` | -9.81 | 重力加速度 y成分 [m/s^2] |
| `gravity_z` | 0.0 | 重力加速度 z成分 [m/s^2] |

### 時間設定

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `dt` | 0.0005 | 時間刻み [s] |
| `t_end` | 2.0 | 終了時刻 [s] |
| `output_interval` | 100 | 出力間隔 [ステップ] |

### 圧力ソルバー

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `solver_type` | 0 | 線形ソルバー (0: CG法, 1: ICCG法) |
| `cg_max_iter` | 10000 | 最大反復回数 |
| `cg_tolerance` | 1.0e-8 | 収束判定閾値 |
| `relaxation_coeff` | 0.2 | 圧力の緩和係数 |

### 自由表面・壁境界

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `surface_threshold` | 0.97 | 自由表面の判定閾値 (n/n0) |

### 計算領域

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `domain_x_min` | 0.0 | 領域 x最小座標 [m] |
| `domain_x_max` | 0.4 | 領域 x最大座標 [m] |
| `domain_y_min` | 0.0 | 領域 y最小座標 [m] |
| `domain_y_max` | 0.3 | 領域 y最大座標 [m] |
| `domain_z_min` | 0.0 | 領域 z最小座標 [m] |
| `domain_z_max` | 0.2 | 領域 z最大座標 [m] |

### 出力

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `output_dir` | output | 出力ディレクトリ |

## 粒子ファイルの形式

`particles.txt` の形式:

```
# コメント行
粒子数
# x y z vx vy vz type
1.0e-02 2.0e-02 3.0e-02 0.0e+00 0.0e+00 0.0e+00 0
...
```

`type`: 0 = 流体粒子, 1 = 壁粒子, 2 = ダミー粒子

## ライセンス

MIT License
