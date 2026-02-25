# MPS法 2次元非圧縮性流体シミュレーション

Moving Particle Semi-implicit (MPS) 法に基づく2次元非圧縮性流体シミュレータ。ダムブレイク問題（水柱崩壊）を対象とし、流体の自由表面流れを粒子法で計算する。

参考文献: S. Koshizuka and Y. Oka, "Moving-Particle Semi-implicit Method for Fragmentation of Incompressible Fluid," Nuclear Science and Engineering, Vol. 123, pp. 421-434, 1996.

## 動作環境

**必須:**
- GCC (または C コンパイラ)
- make
- libm (数学ライブラリ)

**可視化 (オプション):**
- Python 3, NumPy, Matplotlib

## ビルドと実行

```bash
# 初期条件の生成
python3 tools/gen_dambreak.py

# ビルド
make

# 実行
./mps_sim examples/dambreak/cal.txt

# または
make run

# クリーン
make clean
```

実行すると `output/` ディレクトリに CSV および VTK ファイルが出力される。

## 可視化

```bash
python3 vis/visualize.py              # 対話的表示
python3 vis/visualize.py --save       # PNG画像として保存
python3 vis/visualize.py --fluid_only # 流体粒子のみ表示
```

## パラメータ変更

全パラメータは `examples/dambreak/params.txt` で設定する。変更後の再ビルドは不要。

### 粒子設定

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `particle_distance` | 0.008 | 粒子間距離 [m] |
| `influence_ratio_lap` | 2.1 | ラプラシアン用影響半径 / 粒子間距離 |
| `influence_ratio_n` | 2.1 | 粒子数密度用影響半径 / 粒子間距離 |
| `max_neighbors` | 512 | 1粒子あたりの最大近傍数 |
| `wall_layers` | 4 | 壁粒子の層数 |
| `dummy_layers` | 4 | ダミー粒子の層数 |

### 物性値

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `density` | 1000.0 | 密度 [kg/m^3] |
| `viscosity` | 0.0 | 動粘性係数 [m^2/s] |
| `gravity_x` | 0.0 | 重力加速度 x成分 [m/s^2] |
| `gravity_y` | -9.81 | 重力加速度 y成分 [m/s^2] |

### 時間設定

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `dt` | 0.0005 | 時間刻み [s] |
| `t_end` | 2.0 | 終了時刻 [s] |
| `output_interval` | 100 | 出力間隔 [ステップ] |

### 圧力ソルバー

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `solver_type` | 1 | 線形ソルバー (0: CG法, 1: ICCG法) |
| `cg_max_iter` | 10000 | 最大反復回数 |
| `cg_tolerance` | 1.0e-8 | 収束判定閾値 |
| `relaxation_coeff` | 0.2 | 圧力の緩和係数 |
| `clamp_negative_pressure` | 1 | 負圧クランプ (0: 無効, 1: 有効) — PPE求解後に流体粒子の P < 0 を 0 に置換する。引張不安定性の抑制に有効。 |

### 自由表面・壁境界

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `surface_threshold` | 0.97 | 自由表面の判定閾値 (n/n0) |
| `restitution_coeff` | 0.2 | 粒子間衝突の反発係数 |
| `collision_distance_ratio` | 0.5 | 衝突判定距離の係数 (col_dist = ratio × l0) |

### 計算領域

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `domain_x_min` | 0.0 | 領域 x最小座標 [m] |
| `domain_x_max` | 0.584 | 領域 x最大座標 [m] |
| `domain_y_min` | 0.0 | 領域 y最小座標 [m] |
| `domain_y_max` | 0.292 | 領域 y最大座標 [m] |

### 出力

| パラメータ | デフォルト値 | 説明 |
|---|---|---|
| `output_dir` | output | 出力ディレクトリ |

## ファイル構成

```
mps_fluid2D/
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
│   └── visualize.py             # 結果の2D可視化スクリプト
├── examples/
│   └── dambreak/                # ダムブレイク問題の入力ファイル
│       ├── cal.txt              # 計算制御ファイル
│       ├── params.txt           # パラメータファイル
│       └── particles.txt        # 粒子初期条件ファイル
├── Makefile
└── README.md
```

## 粒子ファイルの形式

`particles.txt` の形式:

```
# コメント行
粒子数
# x y vx vy type
1.0e-02 2.0e-02 0.0e+00 0.0e+00 0
...
```

`type`: 0 = 流体粒子, 1 = 壁粒子, 3 = ダミー粒子

## ライセンス

MIT License
