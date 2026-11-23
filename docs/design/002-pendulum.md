# 002: 単振り子

- 状態: 検証済み
- 作成日: 2026-11-23
- 教科書: [第1-1章 単振り子](../textbook/01-mechanics/01-pendulum.md)
- 実装（層ごと、構成は [003](003-layered-architecture.md)）:
  - 物理: [domain/mechanics/pendulum.py](../../src/physlab/domain/mechanics/pendulum.py)
  - ユースケース: [application/mechanics/pendulum.py](../../src/physlab/application/mechanics/pendulum.py)
  - 描画: [infrastructure/rendering/mechanics/pendulum.py](../../src/physlab/infrastructure/rendering/mechanics/pendulum.py)
  - コマンド: [presentation/commands/mechanics/pendulum.py](../../src/physlab/presentation/commands/mechanics/pendulum.py)
- 実行: `uv run physlab run pendulum`

## 1. 目的

- 振れ角が大きいと**小角近似（sin θ ≈ θ）が破綻し、周期が伸びる**様子を動画で見る。
- 数値積分法（オイラー / RK4 / ベルレ）で**エネルギー保存の質がどう違うか**を見る。

## 2. 物理モデル

質量 m の質点、長さ L の軽い剛体棒、重力加速度 g、摩擦なし。鉛直下向きからの角度を θ とする。

```
d²θ/dt² = -(g/L) sin θ            （厳密な運動方程式）
d²θ/dt² = -(g/L) θ                （小角近似）→ θ(t) = θ0 cos(ω0 t),  ω0 = √(g/L)
```

- 力学的エネルギー（単位質量あたり）: `E/m = ½ L² ω² + g L (1 − cos θ)`（最下点を基準）
- 厳密な周期（初期角速度 0）: `T = (4/ω0) K(sin²(θ0/2))`。K は第 1 種完全楕円積分。
  θ0 → 0 で `T → 2π/ω0`。θ0 = 120° では `T ≈ 1.37 × 2π/ω0`。

## 3. 数値手法

| 手法 | 次数 | 用途 |
|---|---|---|
| RK4 | 4 | 動画の本計算（刻み 1 ms） |
| 前進オイラー | 1 | エネルギー比較の対照 |
| 速度ベルレ | 2（シンプレクティック） | エネルギー比較 |

- 動画: dt = 1 ms の RK4。フレームは時刻 k/fps に最も近いステップを使い、再生が実時間と一致する。
- エネルギー比較: 3 手法とも dt = 0.01 s、100 s。

## 4. パラメータと入出力

| 引数 | 既定値 | 意味 |
|---|---|---|
| `--theta0-deg` | 120 | 初期角 [deg]（静かに放す） |
| `--length` | 1.0 | 長さ [m] |
| `--duration` | 10 | 動画の長さ [s] |
| `--fps` | 30 | フレームレート |
| `--energy-duration` | 100 | エネルギー比較の長さ [s] |
| `--energy-dt` | 0.01 | エネルギー比較の刻み [s] |

出力: `outputs/mechanics/pendulum.mp4`, `outputs/mechanics/pendulum_energy_drift.png`

## 5. 可視化の設計

- 動画（左）: 厳密解（青）と小角近似（橙）の 2 つの振り子を重ねて描く。
  90° を超えると錘が支点より上に来るので、表示範囲は初期角から自動で決める。
- 動画（右）: θ(t) の全体を薄く先に描き、現在位置の点と軌跡の太線を動かす。
  タイトルに厳密周期と近似周期を併記する。
- 静止画: 積分法ごとの |E − E0|/E0 を対数軸で比較する（符号は捨てる）。

## 6. 検証方法

層ごとに分けて検証する。

- domain（式・極限・特殊値。時間積分は使わない）:
  - 小角近似の周期が 2π√(L/g)、θ0 → 0 で厳密な周期がそれに一致、120° で T/T0 = 1.3729
  - 周期が振幅とともに単調増加、級数 1 + θ0²/16 との一致（20°）
  - エネルギー式（放した位置と最下点）、運動方程式の特殊点（最下点・水平・真上）
- application（時間積分した結果を照合）:
  - 数値解のゼロ交差から求めた周期が楕円積分の厳密解と相対誤差 1e-4 以内で一致（θ0 = 10, 60, 120, 170°）
  - 小角極限: θ0 = 1° で軌道が cos に一致（最大誤差 1e-4 rad 未満）。120° では大きくずれる
  - エネルギー: RK4・ベルレで 100 s 後も相対変化が 1% 未満。手法間の順序（RK4 < ベルレ < 1e-3 < オイラー）
- numerics: 調和振動子で収束次数（オイラー 1、RK4 4、ベルレ 2）と、オイラーの (1 + dt²) 倍則（`tests/numerics/test_integrators.py`）
- infrastructure / presentation: 煙テスト（ファイルができる）、CLI の入出力
- 全体: 30 件、すべて成功（`uv run pytest`）。移行（003）の前後で実行結果の数値が一致することも確認済み

## 7. 結果・知見

θ0 = 120°, L = 1 m, g = 9.81 m/s²:

- 周期は厳密 2.7541 s、小角近似 2.0061 s。**約 37% 伸びる**。動画では 2 秒付近で橙は一周して
  振り出し位置に戻るが、青はまだ戻れない。
- エネルギー（dt = 0.01 s, 100 s）: オイラー +321%（増え続ける）、RK4 −7.1e-8、ベルレ 振幅 1.8e-4 の有界振動。

詳細は [knowledge/numerical-integration.md](../knowledge/numerical-integration.md) と
[knowledge/visualization.md](../knowledge/visualization.md)。

## 8. 未解決・次にやること

- 位相空間（θ–ω 平面）の軌道図を足すと、シンプレクティック法の性質が見やすくなりそう。
- 長時間（10⁵ s 級）での RK4 とベルレの逆転は未検証（[知見](../knowledge/numerical-integration.md)参照）。
- 減衰・外力を加えた版は「調和振動子（減衰・強制）」の題材で扱う。
