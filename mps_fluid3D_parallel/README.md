# MPS 3D Fluid Simulation (OpenMP並列化版)

MPS (Moving Particle Semi-implicit) 法による3次元非圧縮性流体シミュレーションプログラムの **OpenMP並列化版**です。
`mps_fluid3D` (逐次版) をベースに、主要な計算ループをOpenMPで並列化しています。

## 概要

MPS法は粒子法の一種で、流体をラグランジュ的に追跡する粒子で離散化し、ナビエ・ストークス方程式を解きます。本プログラムでは以下の機能を実装しています。

- 3次元非圧縮性流体の数値シミュレーション
- 半陰的圧力解法 (CG法 / ICCG法を選択可能)
- 自由表面の判定 (粒子数密度 / 近傍粒子数に基づく)
- 壁境界処理 (壁粒子による反発力モデル)
- 表面張力モデル (Colagrossi型ポテンシャル)
- セルリンクリストによる高速近傍粒子探索
- CSV形式 / VTK形式での結果出力
- **OpenMPによるマルチスレッド並列化**

### 圧力計算の対象粒子

| 粒子種別 | 圧力ポアソン方程式 | 圧力勾配による速度・位置更新 |
|---|---|---|
| 流体粒子（内部） | 求解 | あり |
| 流体粒子（自由表面） | 0固定（ディリクレ境界条件） | あり |
| 壁粒子 | 求解 | なし |
| ダミー粒子 | 0固定・圧力計算に非寄与 | なし |

## 並列化の詳細

### 並列化した箇所

| ファイル | 対象処理 | 並列化方法 |
|---|---|---|
| `neighbor_search.c` | 近傍粒子探索 フェーズ2 | `parallel for schedule(dynamic,32)` |
| `operators.c` | 粒子数密度・粘性・圧力勾配・表面張力・衝突 | `parallel for schedule(static/dynamic)` |
| `pressure_solver.c` | CSR行列構築・行列ベクトル積・CG/ICCGのsaxpy | `parallel for schedule(static)` |
| `simulation.c` | 加速度初期化・速度・位置更新 | `parallel for schedule(static)` |
| `boundary.c` | 壁境界・領域外粒子処理 | `parallel for schedule(static)` |

### 逐次のまま維持した箇所

| 処理 | 理由 |
|---|---|
| セルリスト構築 (フェーズ1) | `head[]` へのprepend操作で競合が発生するため |
| 未知数マッピング `eq_idx[]` | prefix-sum のカウントが逐次依存のため |
| IC(0) 分解 | 前の行の計算結果に依存するため |
| 前進/後退代入 | データ依存があるため |
| ドット積 `dot()` | **下記の数値安定性の注記を参照** |

### 数値安定性に関する注記 — `dot()` を逐次にする理由

OpenMP の `reduction(+:s)` はスレッドごとの部分和を後から合計するため、浮動小数点演算の加算順序が逐次版と異なります。CG/ICCG ソルバーでは `alpha = (r·r) / (p·Ap)` などの係数がドット積から計算されるため、この微小な差異が数百ステップにわたって蓄積し、最終的に圧力解が発散するケースが確認されました（テスト: 4スレッド、t≈0.25s で全流体粒子が領域外に飛び出す）。

対策として `dot()` のみ逐次実装とし、逐次版とビット同一の収束挙動を保証しています。行列ベクトル積 `csr_mat_vec()` は各行の内部sumが逐次のまま並列化されるためビット同一となり、saxpy は各要素が独立計算なのでビット同一です。

### 性能計測結果 (dambreakケース、Apple M系10コア)

| 設定 | 実時間 | 速度比 |
|---|---|---|
| 逐次版 (`mps_fluid3D`) | 2:45 | 1.0× |
| 2スレッド | 2:28 | 1.12× |
| **4スレッド（推奨）** | **1:57** | **1.42×** |
| 8スレッド | 2:05 | 1.33× |

4スレッドがベストです。粒子数が多い問題（l0を小さくするなど）ほど並列化効率が高まります。

## ディレクトリ構成

```
mps_fluid3D_parallel/
├── src/                         # ソースコード (.c)
│   ├── main.c                   # メインルーチン・粒子ファイル読み込み
│   ├── simulation.c             # シミュレーションのメインループ (並列化)
│   ├── particle.c               # 粒子システムの管理
│   ├── neighbor_search.c        # 近傍粒子探索 (フェーズ2並列化)
│   ├── kernel.c                 # MPS法のカーネル関数
│   ├── operators.c              # 勾配・ラプラシアン等の微分演算子 (並列化)
│   ├── pressure_solver.c        # 圧力ポアソン方程式ソルバー (部分並列化)
│   ├── boundary.c               # 壁境界処理 (並列化)
│   ├── sim_config.c             # 設定の読み込み・管理
│   └── io.c                     # ファイル出力 (CSV / VTK)
├── include/                     # ヘッダファイル (.h)
├── .clangd                      # clangd用OpenMPインクルードパス設定
├── obj/                         # オブジェクトファイル (ビルド時に生成)
├── output/                      # シミュレーション結果の出力先
├── tools/
│   ├── gen_dambreak.py          # ダムブレイク問題の初期条件生成
│   └── gen_droplet.py           # 液滴振動問題の初期条件生成
├── vis/
│   └── visualize.py             # 結果の3D可視化スクリプト
├── examples/
│   ├── dambreak/                # ダムブレイク問題の入力ファイル
│   │   ├── cal.txt
│   │   ├── params.txt
│   │   └── particles.txt
│   └── droplet/                 # 液滴振動問題の入力ファイル
│       ├── cal.txt
│       ├── params.txt
│       └── particles.txt
└── Makefile
```

## 必要環境

- GCC (C99以降) または Apple Clang
- Make
- **OpenMP ライブラリ**
  - Linux (GCC): `-fopenmp` のみで可 (追加インストール不要)
  - macOS (Apple Clang): `brew install libomp` が必要
- Python 3 (初期条件生成・可視化用)
  - NumPy, Matplotlib

## ビルド

```bash
make
```

Makefile が macOS / Linux を自動判定して適切なフラグを設定します。実行ファイル `mps_sim` が生成されます。

```bash
make clean   # ビルド成果物の削除
```

## 使い方

### 1. 初期条件の生成

```bash
python3 tools/gen_dambreak.py
```

デフォルトでは `examples/dambreak/` に以下の3ファイルが生成されます。

| ファイル | 内容 |
|---|---|
| `cal.txt` | 計算制御ファイル (粒子ファイル・パラメータファイルのパスを指定) |
| `params.txt` | シミュレーションパラメータ |
| `particles.txt` | 粒子の初期位置・速度・種別 |

主要なオプション:

```bash
python3 tools/gen_dambreak.py \
  --l0 0.008 \          # 粒子間距離 [m]（小さいほど精度向上、粒子数増加）
  --domain_x 0.584 \    # 領域 x幅 [m]
  --domain_y 0.292 \    # 領域 y高さ [m]
  --domain_z 0.024 \    # 領域 z奥行 [m]
  --dt 5.0e-4 \         # 時間刻み [s]
  --t_end 2.0           # 終了時刻 [s]
```

### 2. シミュレーションの実行

スレッド数を `OMP_NUM_THREADS` 環境変数で指定します（省略時はシステム最大コア数）。

```bash
# スレッド数を指定して実行（推奨）
OMP_NUM_THREADS=4 ./mps_sim examples/dambreak/cal.txt

# make run でも実行可能（OMP_NUM_THREADSが設定されていればそれを使用）
OMP_NUM_THREADS=4 make run
```

起動時に使用スレッド数が表示されます：

```
OpenMP: using 4 thread(s)
Starting simulation (3D, OpenMP): 4000 steps, dt = 5.00e-04
Step    100 / 4000  (t = 0.0500 s)  fluid particles: 1296
...
```

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

VTKファイルは [ParaView](https://www.paraview.org/) でも開けます。

## パラメータ

`params.txt` で設定可能なパラメータ:

### 粒子設定

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `particle_distance` | 0.008 | 粒子間距離 [m] |
| `influence_ratio_lap` | 2.1 | ラプラシアン用影響半径 / 粒子間距離 |
| `influence_ratio_n` | 2.1 | 粒子数密度用影響半径 / 粒子間距離 |
| `max_neighbors` | 512 | 1粒子あたりの最大近傍数 |
| `wall_layers` | 2 | 壁粒子の層数 |
| `dummy_layers` | 2 | ダミー粒子の層数 |

### 物性値

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `density` | 1000.0 | 密度 [kg/m³] |
| `viscosity` | 1e-06 | 動粘性係数 [m²/s] |
| `gravity_x` | 0.0 | 重力加速度 x成分 [m/s²] |
| `gravity_y` | -9.81 | 重力加速度 y成分 [m/s²] |
| `gravity_z` | 0.0 | 重力加速度 z成分 [m/s²] |

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
| `relaxation_coeff` | 0.2 | 圧力の緩和係数 (`ppe_type=0` のみ) |
| `clamp_negative_pressure` | 0 | 負圧クランプ (0: 無効, 1: 有効) |
| `ppe_type` | 0 | PPE定式化 (0: 既存密度型, 1: Natsui弱圧縮型) |
| `c_ppe` | 1.01 | Natsui型PPEの対角係数 c (`ppe_type=1` のみ) |
| `gamma_ppe` | 0.01 | Natsui型PPEの密度補正重み γ (`ppe_type=1` のみ) |
| `hs_mode` | 0 | 高次ソース項 (0: 無効, 1: 有効) |

### 自由表面判定

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `surface_threshold` | 0.97 | 判定閾値 n/n0 (`surface_detection_method=0` のみ) |
| `surface_detection_method` | 0 | 判定方法 (0: 粒子数密度, 1: 近傍粒子数 Natsui法) |
| `surface_count_threshold` | 0.85 | 近傍粒子数比の閾値 Ni/N0 (`surface_detection_method=1` のみ) |

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
| `domain_x_max` | 0.584 | 領域 x最大座標 [m] |
| `domain_y_min` | 0.0 | 領域 y最小座標 [m] |
| `domain_y_max` | 0.292 | 領域 y最大座標 [m] |
| `domain_z_min` | 0.0 | 領域 z最小座標 [m] |
| `domain_z_max` | 0.024 | 領域 z最大座標 [m] |

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
