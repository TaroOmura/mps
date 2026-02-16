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
make            # ビルド
make run        # 実行
make clean      # クリーン
```

実行すると `output/` ディレクトリに CSV および VTK ファイルが出力される。

## 可視化

```bash
python vis/visualize.py          # 対話的表示
python vis/visualize.py --save   # PNG画像として保存
```

## パラメータ変更

全パラメータは `include/config.h` に定義されている。変更後は `make clean && make` で再ビルドが必要。

| パラメータ | 説明 | デフォルト値 |
|---|---|---|
| `PARTICLE_DISTANCE` | 粒子間距離 | 0.025 m |
| `DT` | 時間刻み | 5.0e-4 s |
| `T_END` | 終了時刻 | 2.0 s |
| `OUTPUT_INTERVAL` | 出力間隔 | 100 ステップ |

## ファイル構成

```
include/    ヘッダファイル (config, particle, simulation, operators, etc.)
src/        ソースファイル (main, simulation, operators, pressure_solver, etc.)
vis/        可視化スクリプト (visualize.py)
output/     シミュレーション出力 (自動生成)
```
