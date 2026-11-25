# 006: 保存法則と力学的相似（ランダウ＝リフシッツ第 II 章の話題）

- 状態: 実装済み・検証済み
- 作成日: 2026-11-26
- 教科書: [第 4 章 保存法則](../textbook/tex/chapters/02-analytical-mechanics/02-conservation-laws.tex)
- 実装（層ごと、構成は [003](003-layered-architecture.md)）:
  - 物理: `domain/mechanics/{cluster,rotating_rod,central,power_law}.py`
  - ユースケース: `application/mechanics/{symmetry,center_of_mass,areal_velocity,similarity}.py`
  - 描画: `infrastructure/rendering/mechanics/{symmetry,center_of_mass,areal_velocity,similarity}.py`
  - コマンド: `presentation/commands/mechanics/{symmetry,center_of_mass,areal_velocity,similarity}.py`
- 実行: `uv run physlab run symmetry` / `center-of-mass` / `areal-velocity` / `similarity`

## 1. 目的と構成の方針

扱う話題は §6 エネルギー、§7 運動量、§8 慣性中心、§9 角運動量、§10 力学的相似。
参考書の本文は写さず、[著作権の方針](../knowledge/copyright.md)に従って構成を自分で決める。
軸は「**対称性を 1 つずつ壊すと、どの量が保存しなくなるか**」を数値実験で見ること。

| 見えれば成功 | コマンド | 図・動画 |
|---|---|---|
| 対称性と保存量の対応（壊した対称性に対応する量だけが変化し、変化の速さが予測と合う） | `symmetry` | `cons-symmetry`、`symmetry_cluster.mp4` |
| 数値解法はどの保存量を守るか（ベルレと RK4） | `symmetry` | `cons-integrators` |
| 時間に依存する拘束があると $\sum p\dot q-L$ は保存しても $T+U$ は保存しない | `symmetry` | `cons-rotating-rod` |
| 投げた物体は回って振動しても、重心は放物線。内部エネルギーと重心まわりの角運動量は保存 | `center-of-mass` | `com-throw`、`com-decomposition`、`com_throw.mp4` |
| 中心力なら面積速度は一定（同じ時間に同じ面積） | `areal-velocity` | `am-areal`、`areal_velocity.mp4` |
| ポテンシャルが同次式なら、運動は相似に写り合う。周期の冪、重ね書き、ビリアル定理 | `similarity` | `sim-periods`、`sim-collapse`、`sim-virial` |

## 2. 物理モデル

### 2.1 ばねでつないだ 3 質点（`cluster`）

$m=(1,2,1.5)$、全ペアを自然長 $\ell=1$、ばね定数 $k$ のばねで結ぶ（2 次元）。
内力は 2 質点の相対位置だけに依り、中心力。外場 $U_{\rm ext}$ を 4 通りにする。

| 外場 | $U_{\rm ext}$ | 壊れる対称性 | 保存する量 |
|---|---|---|---|
| `none` | 0 | なし | $E,\ P_x,\ P_y,\ L$ |
| `gravity` | $\sum m_agy_a$ | $y$ の並進、回転 | $E,\ P_x$ |
| `trap` | $\tfrac12\kappa\sum m_a r_a^2$ | 並進 | $E,\ L$ |
| `pulsing_trap` | $\tfrac12\kappa(t)\sum m_ar_a^2$、$\kappa=\kappa_0(1+\beta\sin\Omega t)$ | 並進、時間の一様性 | $L$ |

保存しない量の変化の速さ（予測）: $\dot{\bm P}=-\sum_a\partial U_{\rm ext}/\partial\bm r_a$、
$\dot L=\sum_a(\bm r_a\times\bm F^{\rm ext}_a)_z$、$\dot E=\partial U_{\rm ext}/\partial t$。
これを時間積分した値と、数値解の $E,\bm P,L$ の変化を比べる。

### 2.2 回転する棒の上のビーズ（`rotating_rod`）

一定の角速度 $\omega$ で回る水平な棒に通した質点。$\mathcal L=\tfrac12m(\dot r^2+\omega^2r^2)$（時刻を陽に含まない）。
$E=p_r\dot r-\mathcal L=\tfrac12m(\dot r^2-\omega^2r^2)$ は保存し、$T=\tfrac12m(\dot r^2+\omega^2r^2)$ は増える。
$r(0)=r_0,\ \dot r(0)=0$ の解は $r=r_0\cosh\omega t$。棒を回し続けるにはモーターが仕事率 $2mr\dot r\,\omega^2$
を与える必要があり、その仕事が $T$ の増加にちょうど等しい：
$T(t)-T(0)=m\omega^2r_0^2\sinh^2\omega t=\int 2mr\dot r\,\omega^2\,dt$（数値で確かめる）。

### 2.3 一様重力に投げた 3 質点（`center_of_mass`）

2.1 の系（$k$ を大きくする）に重心の初速と回転を与える。$g=9.81$。
期待: $\bm R(t)=\bm R_0+\bm V_0t+\tfrac12\bm gt^2$、内部エネルギー $E_{\rm int}=T_{\rm rel}+U_{\rm spring}$ と、
重心まわりの角運動量 $L_{\rm rel}$ は保存する。

### 2.4 中心力（`central`）

$U=-\alpha/r^n$（$0<n<2$）の 2 次元運動。$n=1.4$ でロゼット状の軌道。
面積速度 $\dot S=L/(2m)$ は一定。速度ベルレ法では 1 ステップの三角形の面積 $\tfrac12|\bm r_k\times\bm r_{k+1}|$ が
$\tfrac{h}{2m}|L|$ に**厳密に**一致する（加速度が $\bm r_k$ に平行なため）。

### 2.5 冪ポテンシャル（`power_law`）

1 次元、$m=1$、$U=|x|^k/k$（$k=1,2,4,6$）。振幅 $A$ から静かに放すときの周期は厳密に

$$T(A)=4\sqrt{k/2}\,A^{1-k/2}\,\frac{\Gamma(1/k)\,\Gamma(1/2)}{k\,\Gamma(1/k+1/2)}$$

力学的相似: $x_A(t)=A\,x_1(tA^{k/2-1})$。ビリアル定理: 時間平均で $2\langle T\rangle=k\langle U\rangle$。
$k=-1$（万有引力）は二体問題（`two_body`）で、$T\propto a^{3/2}$ と $2\langle T\rangle=-\langle U\rangle$ を確かめる。

## 3. 数値手法

| 処理 | 手法 | 刻み |
|---|---|---|
| 3 質点（対称性の実験、投げ上げ） | 速度ベルレ（比較として RK4） | $h=1$ ms（比較は 10 ms） |
| 中心力 | 速度ベルレ | $h=0.5$ ms |
| 冪ポテンシャルの周期 | 速度ベルレ、ゼロ交差の線形補間 | $T_{\rm exact}/4000$ |

## 4. パラメータと入出力

| コマンド | 主な引数 | 出力 |
|---|---|---|
| `symmetry` | `--duration 10`、`--no-video` | `outputs/mechanics/symmetry_cluster.mp4`、図 3 枚 |
| `center-of-mass` | `--no-video` | `outputs/mechanics/com_throw.mp4`、図 2 枚 |
| `areal-velocity` | `--exponent 1.4`、`--no-video` | `outputs/mechanics/areal_velocity.mp4`、図 1 枚 |
| `similarity` | なし | 図 3 枚 |

- `center-of-mass` の動画は、飛行が約 1.6 s と短いので 1/4 倍速で再生し、倍率を画面に表示する
  （他の動画は実時間）。
- 教科書の図（`figures()`）は引数を取らず、常に上の既定の条件で作る。

## 5. 検証方法（必須）

| 層 | 検証 | 期待 |
|---|---|---|
| domain | 3 質点の力 $=-\nabla U$（数値微分）、内力の和と内力のモーメントが 0 | $10^{-6}$ / 丸め |
| domain | 外場ごとの $\dot E,\dot{\bm P},\dot L$ の式が、運動方程式からの連鎖律と一致 | $10^{-9}$ |
| domain | 冪ポテンシャルの周期の式: $k=2$ で $2\pi$、$k=1$ で $4\sqrt{2A}$ | $10^{-12}$ |
| domain | 回転棒の解が運動方程式を満たし $h$ が一定 | $10^{-9}$ |
| application | 外場ごとに保存する量は丸め／積分誤差の範囲で一定、他は予測どおり変化 | 予測との差 $10^{-5}$ |
| application | ベルレは $\bm P,L$ を丸め誤差で保存。RK4 は $\bm P$ は保存、$L$ はずれる | $10^{-12}$ / 有意な差 |
| application | 重心が放物線、$E_{\rm int}$ と $L_{\rm rel}$ が保存、$M\bm R-\bm Pt$（自由な系）が一定 | $10^{-10}$ / $10^{-5}$ |
| application | 同じ時間の扇形の面積が等しい。弧の長さは大きく異なる | 相対 $10^{-12}$ |
| application | 周期が公式と一致、傾きが $1-k/2$、相似変換した軌道が重なる、$2\langle T\rangle/\langle U\rangle=k$ | $10^{-4}$ |
| 目視 | PDF の該当ページと動画のフレーム | 切れ・重なり・文字化けなし |

## 6. 結果・知見

実装後に [knowledge/conservation-laws.md](../knowledge/conservation-laws.md) に実測値を書く。

## 7. 未解決・次にやること

- 3 次元の角運動量（ベクトルとしての保存）は扱わない（2 次元の $z$ 成分のみ）。
- ネーターの定理の場の理論版には触れない。
