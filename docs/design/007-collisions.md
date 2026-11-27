# 007: 粒子の衝突と散乱（ランダウ＝リフシッツ第 IV 章の話題）

- 状態: 実装済み・検証済み
- 作成日: 2026-11-27
- 教科書: [第 5 章 粒子の衝突と散乱](../textbook/tex/chapters/02-analytical-mechanics/03-collisions.tex)
- 実装（層ごと、構成は [003](003-layered-architecture.md)）:
  - 物理: `domain/mechanics/{decay,collision,potentials,scattering,alpha_gold}.py`
  - 数値: `numerics/{quadrature,roots}.py`
  - ユースケース: `application/mechanics/{decay,elastic_collision,deflection,scattering,rutherford,small_angle}.py`
  - 描画: `infrastructure/rendering/mechanics/{decay,elastic_collision,scattering,rutherford,small_angle}.py`
  - コマンド: `presentation/commands/mechanics/{decay,elastic_collision,scattering,rutherford,small_angle}.py`
- 実行: `uv run physlab run decay` / `elastic-collision` / `scattering` / `rutherford` / `small-angle`

## 1. 目的と構成の方針

扱う話題は §16 粒子の崩壊、§17 弾性衝突、§18 散乱、§19 ラザフォードの公式、§20 微小角散乱。
参考書の本文は写さず、[著作権の方針](../knowledge/copyright.md)に従って構成を自分で決める。
軸は「**力を知らなくても分かること**（保存則）と、**力を知らないと分からないこと**（散乱角）」。

第 III 章（運動方程式の積分）の章はまだ無いので、中心力の軌道の積分は、この章の中で
エネルギーと角運動量の保存から必要な分だけ導く。

| 見えれば成功 | コマンド | 図・動画 |
|---|---|---|
| 崩壊した破片の実験室系の速度は 1 つの球面上にある。親が速いと円すいの中だけに飛ぶ | `decay` | `dec-velocity`、`dec-distributions`、`decay_cone.mp4` |
| 力の形がまったく違っても、衝突後の運動量は同じ円の上。違うのは円の上の位置だけ | `elastic-collision` | `col-momentum`、`col-trajectories`、`col-mass-ratio`、`elastic_collision.mp4` |
| 衝突パラメータ → 散乱角 → 断面積。乱数で打ち込んで数えた分布が理論と合う。虹散乱 | `scattering` | `sca-deflection`、`sca-cross-section`、`sca-lab-frame`、`scattering_beam.mp4` |
| クーロン力の断面積は $\sin^{-4}(\chi/2)$。電荷を広げると大角散乱が消える | `rutherford` | `rut-orbits`、`rut-cross-section`、`rut-charge-size`、`rutherford.mp4` |
| 直線の経路に沿った力積で散乱角が求まる。誤差の減り方 | `small-angle` | `sma-impulse`、`sma-error` |

## 2. 物理モデル

記号: 換算質量 $\mu=m_1m_2/(m_1+m_2)$、重心系の散乱角 $\chi$、実験室系の角 $\theta_1,\theta_2$、
衝突パラメータ $\rho$、無限遠の速さ $v_\infty$、衝突のエネルギー $E=\tfrac12\mu v_\infty^2$。
単位系は無次元（計算では長さ・質量・速さの基準を 1 にする）。

### 2.1 崩壊（`decay`）

内部エネルギー $\varepsilon$ が解放されて 2 つに分かれる。親の静止系で $p_0=\sqrt{2\mu\varepsilon}$、
$v_{10}=p_0/m_1$。親の速度 $\bm V$ の系では $\bm v=\bm V+v_0\bm n$。$k=V/v_0$ とすると

$$\tan\theta=\frac{\sin\theta_0}{\cos\theta_0+k},\qquad \theta_{\max}=\arcsin(1/k)\ (k>1)$$

親の静止系で等方的なとき、実験室系の角の確率密度は
$k<1$: $\tfrac12\sin\theta\,[2k\cos\theta+(1+k^2\cos2\theta)/\sqrt{1-k^2\sin^2\theta}]$、
$k\ge1$: $\sin\theta\,(1+k^2\cos2\theta)/\sqrt{1-k^2\sin^2\theta}$（$\theta<\theta_{\max}$）。
運動エネルギーは $\tfrac12m(V-v_0)^2$ から $\tfrac12m(V+v_0)^2$ まで一様に分布する。

### 2.2 弾性衝突（`collision`）

重心系では相対速度の向きが $\chi$ だけ回る。静止した標的（$\bm v_2=0$）では

$$\tan\theta_1=\frac{m_2\sin\chi}{m_1+m_2\cos\chi},\quad\theta_2=\frac{\pi-\chi}{2},\quad
\frac{E_2'}{E_1}=\frac{4m_1m_2}{(m_1+m_2)^2}\sin^2\frac\chi2$$

衝突後の $\bm p_2'$ は、中心 $(\mu v,0)$・半径 $\mu v$ の円の上にある。
数値実験では、2 粒子の運動方程式を、形の違う短距離の斥力 4 種類で解く
（柔らかい円板 $\tfrac12k(a-r)^2$、ガウス型 $U_0e^{-r^2/a^2}$、湯川型 $\alpha e^{-r/\lambda}/r$、
急な芯 $\epsilon(a/r)^{12}$）。

### 2.3 散乱（`potentials`、`scattering`）

固定した中心力場（相対運動に置き換えた形）で、$u=1/r$ とすると

$$\varphi_0=\int_0^{u_{\max}}\frac{\rho\,\dd u}{\sqrt{F(u)}},\qquad
F(u)=1-\rho^2u^2-\frac{U(1/u)}{E},\qquad \Theta=\pi-2\varphi_0$$

$u_{\max}$ は $F$ の最初の零点（最接近）。$\Theta$ は符号つきの偏向角（斥力で正）、$\chi=|\Theta|$。
断面積は $\dd\sigma/\dd\Omega=\sum_{\text{枝}}(\rho/\sin\chi)\,|\dd\rho/\dd\chi|$。
剛体球: $\rho=a\cos(\chi/2)$、$\dd\sigma/\dd\Omega=a^2/4$。

**仮定**: $F$ の正の零点が 1 つ（斥力で単調な場合と、クーロン引力）。軌道が中心のまわりを回り続ける
引力（orbiting）は扱わない。

例として使うポテンシャル: クーロン $\alpha/r$（$\alpha>0$ が斥力）、逆 2 乗 $\alpha/r^2$、
ガウス型の丘、一様に帯電した球の外と中。

### 2.4 ラザフォード

$\rho=(d/2)\cot(\chi/2)$、$d=|\alpha|/E$（斥力の正面衝突の最接近距離）、
$\dd\sigma/\dd\Omega=(d/4)^2\sin^{-4}(\chi/2)$、最接近 $r_{\min}=(d/2)(1+1/\sin(\chi/2))$（斥力）。
静止標的への移行エネルギー: $\dd\sigma=(2\pi\alpha^2/m_2v_\infty^2)\,\dd\epsilon/\epsilon^2$。

一様に帯電した球（半径 $R$）: 外は $\alpha/r$、中は $\alpha(3-r^2/R^2)/2R$。
$R<d$ なら点電荷と同じ。$d<R$ では、$\chi<\chi_g=2\arcsin\dfrac{1}{2R/d-1}$ が点電荷と同じで、
それより大きい角で外れる。$R>1.5d$ なら中心を通り抜けられ、散乱角に最大値がある。

### 2.5 微小角散乱

直線 $y=\rho$ に沿った力積から

$$\theta=-\frac{2\rho}{\mu v_\infty^2}\int_\rho^\infty\frac{\dd U}{\dd r}\frac{\dd r}{\sqrt{r^2-\rho^2}},\qquad
U=\frac{\alpha}{r^n}\ \Rightarrow\ \theta=\frac{2\sqrt\pi\,\alpha\,\Gamma\!\left(\frac{n+1}2\right)}
{\mu v_\infty^2\rho^n\,\Gamma\!\left(\frac n2\right)}$$

このとき $\dd\sigma/\dd\Omega\approx(1/n)A^{2/n}\theta^{-2-2/n}$（$\theta=A\rho^{-n}$）。

## 3. 数値手法

| 処理 | 手法 | 刻み・点数 |
|---|---|---|
| 2 体衝突の軌道 | 速度ベルレ（複数の衝突を 1 本の状態ベクトルにまとめて同時に解く） | $h=2\times10^{-3}$ |
| 固定力場の軌道 | 速度ベルレ（粒子の束を同時に解く） | $h=10^{-3}$ |
| 最接近 $u_{\max}$ | 2 分法（配列でまとめて） | 反復 100 回 |
| 偏向角の積分 | ガウス・ルジャンドル求積。区間を 2 つに分け、それぞれ上端に節点を集める置換 $u=b-(b-a)s^2$ を使う（端の $1/\sqrt{\ }$ を消し、球の表面をかすめる軌道の準特異性にも効く） | 120 点 × 2（数える実験は 64 点 × 2） |
| 力積の積分 | $r=\rho/\cos u$ の置換 + ガウス・ルジャンドル | 120 点 |
| 数える実験 | 衝突パラメータを円板内で一様に乱数で選ぶ（$\rho=\rho_{\max}\sqrt{\xi}$）。乱数の種は固定 | $2\times10^5$ 個 |

- 引力でほぼ正面から入る軌道は求積の精度が落ちる（使わない範囲）。
- 軌道の初期条件は、有限の距離 $R_0$ から始めても無限遠の $\rho$ に正しく対応するように、
  入射側の偏角 $\int_0^{1/R_0}\rho\,\dd u/\sqrt F$ を使って置く。終状態の向きも同じ積分で無限遠まで延ばす。
- ポテンシャルに折れ目がある場合（帯電球の表面）は、積分区間をそこで分ける。
- 偏向角は 2 通り（軌道の時間積分と、軌道の式の求積）で計算して、互いに照合する。

## 4. パラメータと入出力

| コマンド | 主な引数 | 出力 |
|---|---|---|
| `decay` | `--samples`、`--no-video` | 図 2 枚、`outputs/mechanics/decay_cone.mp4` |
| `elastic-collision` | `--m2 2.0`、`--impact 0.6`、`--no-video` | 図 3 枚、`elastic_collision.mp4` |
| `scattering` | `--energy 1.5`（$E/U_0$）、`--no-video` | 図 3 枚、`scattering_beam.mp4` |
| `rutherford` | `--ball-radius 1.4`（$R/d$）、`--no-video` | 図 3 枚、`rutherford.mp4` |
| `small-angle` | なし | 図 2 枚 |

引数は画面に出す数値に効く（`elastic-collision` の `--m2`・`--impact` と `scattering` の `--energy` は動画にも）。
教科書の図（`figures()`）は引数を取らず、常に既定の条件で作る。

## 5. 可視化の設計

- 速度の図: 親の速さ $\bm V$ を中心とする半径 $v_0$ の円。$k<1$ と $k>1$ を並べる。
- 運動量の図: 円の上に、4 種類の力 × いくつもの $\rho$ の数値解の終状態を点で重ねる。
- 散乱の動画: 左に力場と粒子の束、右の偏向角の曲線の上に、力の中心を通り過ぎた粒子の $(\rho,\chi)$ が現れる
  （一様な $\rho$ の束は面積に比例した分布ではないので、ヒストグラムにはしない）。
- ラザフォードの断面積: 角度が大きいほど桁で小さくなるので縦軸は対数。統計の誤差棒をつける。

## 6. 検証方法（必須）

| 層 | 検証 | 期待 |
|---|---|---|
| numerics | 求積: 多項式、端点に $1/\sqrt{\ }$ の特異性がある積分の厳密値 / 2 分法: 既知の根 | $10^{-12}$ |
| domain | 崩壊の角度分布の積分が 1（台形公式）。$V=0$ で $\tfrac12\sin\theta$ | $2\times10^{-4}$ |
| domain | 弾性衝突の終速度が運動量とエネルギーを保存。等質量で $\theta_1+\theta_2=90^\circ$。$m_1>m_2$ の最大角 | $10^{-12}$ |
| domain | ポテンシャルの $\dd U/\dd r$ が数値微分と一致。帯電球は表面で値と傾きが連続 | $10^{-6}$ |
| domain | ラザフォード・剛体球・逆 2 乗の式どうしの整合（$\dd\sigma/\dd\Omega=(\rho/\sin\chi)\lvert\dd\rho/\dd\chi\rvert$） | $10^{-6}$ |
| application | 崩壊の乱数実験が運動量・エネルギーを保存し、分布が理論と統計誤差の範囲で一致 | 4σ |
| application | 4 種類の力の衝突の終運動量が円の上。$\theta_1,\theta_2$ が $\chi$ からの式と一致 | $10^{-4}$ |
| application | 求積の偏向角が解析解（クーロン引力・斥力、逆 2 乗）と一致 | $10^{-8}$ |
| application | 求積の偏向角が軌道の時間積分と一致（ガウス型、帯電球、クーロン） | $10^{-5}$ rad |
| application | 数えた断面積が剛体球 $a^2/4$、ラザフォードの式と統計誤差の範囲で一致 | 4σ |
| application | 帯電球: $R<d$ は点電荷と同じ、$d<R$ は $\chi_g$ 以下が同じ、$R>1.5d$ は最大角あり | 一致 |
| application | 力積近似の相対誤差の冪: クーロンは $\theta^2$、逆 2 乗やガウス型は $\theta^1$ | 傾き ±0.05 |
| 目視 | PDF の該当ページと動画のフレーム | 切れ・重なり・文字化けなし |

## 7. 結果・知見

実測値とハマった点は [knowledge/collisions.md](../knowledge/collisions.md)。
主な値: 衝突後の運動量の円からのずれ $2\times10^{-6}$ 以下、求積と時間積分の偏向角の差 $6\times10^{-8}$ rad、
数える実験のずれは区間ごとの最大で 2.0〜3.4σ。

## 8. 未解決・次にやること

- 中心に回り込む引力（orbiting、栄光）は扱わない。
- 非弾性衝突、相対論的な運動学、量子力学の散乱は扱わない。
- 標的が動けるラザフォード散乱の軌道の動画は作らない（重心系への換算で扱う）。
