# MPS 3D Fluid Simulation

MPS (Moving Particle Semi-implicit) 法による3次元非圧縮性流体シミュレーションプログラムです。

## 概要

MPS法は粒子法の一種で、流体をラグランジュ的に追跡する粒子で離散化し、ナビエ・ストークス方程式を解きます。本プログラムでは以下の機能を実装しています。

- 3次元非圧縮性流体の数値シミュレーション
- 半陰的圧力解法 (CG法 / ICCG法を選択可能)
- 自由表面の判定 (粒子数密度 / 近傍粒子数に基づく)
- 壁境界処理 (壁粒子による反発力モデル)
- 表面張力モデル (Colagrossi型ポテンシャル): 液液 (C_LL) および固液濡れモデル (C_SL, Young-Dupré式)
- 近傍粒子探索
- CSV形式 / VTK形式での結果出力

### 圧力計算の対象粒子

| 粒子種別 | 圧力ポアソン方程式 | 圧力勾配による速度・位置更新 |
|---|---|---|
| 流体粒子（内部） | 求解 | あり |
| 流体粒子（自由表面） | 0固定（ディリクレ境界条件） | あり |
| 壁粒子 | 求解 | なし |
| ダミー粒子 | 0固定・圧力計算に非寄与 | なし |

- 壁粒子は圧力ポアソン方程式の未知数として含めることで、流体と壁の境界での圧力を適切に計算する。ただし、圧力勾配による加速度・速度・位置の修正は行わない。
- ダミー粒子は壁粒子の粒子数密度 `n_i` を正しく計算するためにのみ使用され、圧力ポアソン方程式には寄与しない（係数行列の対角・非対角成分いずれにも加算されない）。

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
│   ├── gen_dambreak.py          # ダムブレイク問題の初期条件生成スクリプト
│   └── gen_droplet.py           # 液滴振動問題の初期条件生成スクリプト
├── vis/
│   └── visualize.py             # 結果の3D可視化スクリプト
├── examples/
│   ├── dambreak/                # ダムブレイク問題の入力ファイル
│   │   ├── cal.txt              # 計算制御ファイル
│   │   ├── params.txt           # パラメータファイル
│   │   └── particles.txt        # 粒子初期条件ファイル
│   └── droplet/                 # 液滴振動問題の入力ファイル (gen_droplet.py で生成)
│       ├── cal.txt
│       ├── params.txt
│       └── particles.txt
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
| `influence_ratio_lap` | 2.1 | ラプラシアン用影響半径 / 粒子間距離 |
| `influence_ratio_n` | 2.1 | 粒子数密度用影響半径 / 粒子間距離 |
| `max_neighbors` | 512 | 1粒子あたりの最大近傍数 |
| `wall_layers` | 2 | 壁粒子の層数 |
| `dummy_layers` | 2 | ダミー粒子の層数 |

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
| `relaxation_coeff` | 0.2 | 圧力の緩和係数 (`ppe_type=0` のみ使用) |
| `clamp_negative_pressure` | 0 | 負圧クランプ (0: 無効, 1: 有効) — PPE求解後に流体粒子の P < 0 を 0 に置換する。引張不安定性の抑制に有効。 |
| `ppe_type` | 0 | PPE定式化 (0: 既存密度型, 1: Natsui弱圧縮型) |
| `c_ppe` | 1.01 | Natsui型PPEの対角係数 c (`ppe_type=1` のみ使用) |
| `gamma_ppe` | 0.01 | Natsui型PPEの密度補正重み γ (`ppe_type=1` のみ使用) |
| `hs_mode` | 0 | 高次ソース項 (0: 無効, 1: 有効) |

### 自由表面判定

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `surface_threshold` | 0.97 | 自由表面の判定閾値 n/n0 (`surface_detection_method=0` のみ使用) |
| `surface_detection_method` | 0 | 判定方法 (0: 粒子数密度, 1: 近傍粒子数 Natsui法) |
| `surface_count_threshold` | 0.85 | 近傍粒子数比の閾値 Ni/N0 (`surface_detection_method=1` のみ使用) |

### 衝突モデル

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `restitution_coeff` | 0.2 | 粒子間衝突の反発係数 e (0: 完全非弾性, 1: 完全弾性) |
| `collision_distance_ratio` | 0.5 | 衝突判定距離の係数 (col_dist = ratio × l0) |

### 表面張力

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `surface_tension_enabled` | 0 | 表面張力の有効化 (0: 無効, 1: 有効) |
| `surface_tension_coeff` | 0.0728 | 液液表面張力係数 σ_LL [N/m] |
| `surface_tension_re_ratio` | 3.2 | 表面張力影響半径の倍率 (re_st = ratio × l0) |
| `wetting_angle_SL` | 90.0 | 静的接触角 θ [deg] — 固液表面張力係数 C_SL = 0.5(1+cosθ)×C_LL (Young-Dupré式) で計算される。θ=90° のとき C_SL=0 (濡れなし)。`surface_tension_enabled=1` のとき流体-壁粒子ペアに適用される。 |

### 圧力勾配・λ

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `cmps_gradient` | 0 | 圧力勾配モデル (0: 標準 P_j-P_min, 1: CMPS対称型 P_i+P_j-P_imin-P_jmin (Khayyer & Gotoh 2008), 2: Oochi対称型 P_i+P_j (Oochi 2010)) |
| `use_analytical_lambda` | 0 | λの計算方法 (0: 初期粒子配置から計算, 1: 解析解) |

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

`type`: 0 = 流体粒子, 1 = 壁粒子, 2 = ゴースト粒子, 3 = ダミー粒子

## ライセンス

MIT License
