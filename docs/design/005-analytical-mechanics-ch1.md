# 005: 解析力学入門（ランダウ＝リフシッツ『力学』第 I 章 §1–§5）

- 状態: 検証済み
- 作成日: 2026-11-23
- 教科書: [第 3 章 運動方程式](../textbook/tex/chapters/02-analytical-mechanics/01-equations-of-motion.tex)
- 実装（層ごと、構成は [003](003-layered-architecture.md)）:
  - 物理: `domain/mechanics/{lagrangian,double_pendulum,two_body,galilean}.py`
  - 数値: `numerics/optimize.py`（最急降下法）
  - ユースケース: `application/mechanics/{least_action,generalized_coordinates,galilean}.py`
  - 描画: `infrastructure/rendering/mechanics/{least_action,generalized_coordinates,galilean}.py`
  - コマンド: `presentation/commands/mechanics/{least_action,generalized_coordinates,galilean}.py`
- 実行: `uv run physlab run least-action` / `generalized-coordinates` / `galilean`

## 1. 目的

ランダウ＝リフシッツ『力学』第 I 章（一般座標、最小作用の原理、ガリレイの相対性原理、
自由な質点のラグランジアン、質点系のラグランジアン）を、**式だけでは想像しにくい部分を図と
小さなシミュレーションで見える化**しながら理解する。特に次の 3 点が「見えれば成功」。

1. **ハミルトンの原理**: 同じ両端を結ぶ無数の経路のうち、真の運動だけが作用の停留点になる。
   作用 $S[q]$ を「経路の空間の上の地形」として見る。作用を下げ続けると経路が真の運動に収束する。
2. 「最小」作用は**短い時間でだけ**最小で、長いと鞍点になる（停留であって最小とは限らない）。
3. ガリレイ変換で運動方程式は形を変えず、ラグランジアンは**全微分だけ**変わる
   （作用は端点だけで決まる定数がずれる）。

## 2. 構成

| 教科書の節 | 内容 | コマンド | 図・動画 |
|---|---|---|---|
| §1 一般座標 | 自由度、一般座標、$(q,\dot q)$ で状態が決まる | `generalized-coordinates` | `gc-same-q`, `gc-double-pendulum-config`（図）、`double_pendulum.mp4` |
| §2 最小作用の原理 | 経路の族と $S(\varepsilon)$、最小か鞍点か、作用の地形、経路の緩和 | `least-action` | `la-path-family`, `la-second-variation`, `la-curvature`, `la-landscape`, `la-relaxation`（図）、`least_action_relaxation.mp4` |
| §3 ガリレイの相対性 | 電車の中の投げ上げ、運動方程式の不変性 | `galilean` | `gal-train`（図） |
| §4 自由な質点 | $L'=L+df/dt$、作用のずれは経路によらない | `galilean` | `gal-action-shift`（図） |
| §5 質点系 | 二体問題、重心系と実験室系、$\mathrm{EL}\Leftrightarrow$ ニュートン | `galilean` | `gal-two-body`（図）、`galilean_two_body.mp4` |

## 3. 物理モデル

### 3.1 一般座標（§1）

- 単振り子（既存 `domain/mechanics/pendulum.py`）: 同じ初期角 $\theta_0$ で初期角速度 $\omega_0$ だけを変える。
  エネルギー保存 $\tfrac12 L^2\omega_0^2 - gL\cos\theta_0 = -gL\cos\theta_{\max}$ から、
  振れの最大角が $\cos\theta_{\max}=\cos\theta_0-L\omega_0^2/(2g)$ と決まる（検証に使う）。
- 二重振り子（$m_1,m_2,l_1,l_2$、一般座標 $\theta_1,\theta_2$、$\Delta=\theta_1-\theta_2$）:

  $$T=\tfrac12(m_1+m_2)l_1^2\dot\theta_1^2+\tfrac12m_2l_2^2\dot\theta_2^2+m_2l_1l_2\dot\theta_1\dot\theta_2\cos\Delta,\quad
  U=-(m_1+m_2)gl_1\cos\theta_1-m_2gl_2\cos\theta_2$$

  オイラー＝ラグランジュ方程式（質量行列の形 $M(\theta)\ddot\theta=\boldsymbol F$）:

  $$\begin{pmatrix}(m_1+m_2)l_1^2 & m_2l_1l_2\cos\Delta\\ m_2l_1l_2\cos\Delta & m_2l_2^2\end{pmatrix}
  \begin{pmatrix}\ddot\theta_1\\ \ddot\theta_2\end{pmatrix}=
  \begin{pmatrix}-m_2l_1l_2\dot\theta_2^2\sin\Delta-(m_1+m_2)gl_1\sin\theta_1\\
  m_2l_1l_2\dot\theta_1^2\sin\Delta-m_2gl_2\sin\theta_2\end{pmatrix}$$

  微小振動の基準振動数は、$M_0=M(0)$、$K=\mathrm{diag}((m_1+m_2)gl_1,\ m_2gl_2)$ として
  $\det(K-\omega^2M_0)=0$ の根。$m_1=m_2,\ l_1=l_2=l$ なら $\omega^2=(2\mp\sqrt2)\,g/l$。

### 3.2 最小作用の原理（§2）

1 次元の質点 $\mathcal L=\tfrac12m\dot q^2-V(q)$。端点 $(t_1,q_1),(t_2,q_2)$、$\tau=t_2-t_1$、$s=t-t_1$。

| 系 | $V(q)$ | 停留経路（EL の解） |
|---|---|---|
| 自由粒子 | $0$ | $q_1+(q_2-q_1)s/\tau$ |
| 一様重力（$q$ は高さ） | $mgq$ | $q_1+(q_2-q_1)s/\tau+\tfrac12g\,s(\tau-s)$ |
| 調和振動子（$\omega=\sqrt{k/m}$） | $\tfrac12m\omega^2q^2$ | $\dfrac{q_1\sin\omega(\tau-s)+q_2\sin\omega s}{\sin\omega\tau}$（$\sin\omega\tau\ne0$） |

- **経路の族**: $q_\varepsilon(t)=q_*(t)+\varepsilon\,\eta_n(t)$、$\eta_n=\sin(n\pi s/\tau)$（両端で 0）。
  3 つの系はいずれも $\mathcal L$ が $(q,\dot q)$ の 2 次式なので、$S$ は $\varepsilon$ の**ちょうど 2 次式**になる:

  $$S(\varepsilon)=S_*+c_n\varepsilon^2,\qquad c_n=\frac{m\tau}{4}\left[\left(\frac{n\pi}{\tau}\right)^2-\omega^2\right]$$

  （自由・重力は $\omega=0$。1 次の項は EL 方程式のため 0。）$c_n>0$ なら最小、$c_n<0$ の方向があれば鞍点。
  $c_1$ が負になるのは $\tau>\pi/\omega$（調和振動子、運動学的焦点）。
- **作用の地形**: 2 次元の族 $q_*+a\eta_1+b\eta_2$ の $S(a,b)$。$c_1,c_2$ の符号で碗か鞍かが決まる。
- **離散作用**: 時刻 $t_k=t_1+k\Delta t$（$k=0..N$、両端固定）で

  $$S_d=\sum_{k=0}^{N-1}\Delta t\left[\frac m2\Big(\frac{q_{k+1}-q_k}{\Delta t}\Big)^2-\frac12\big(V(q_k)+V(q_{k+1})\big)\right]$$

  内部の節点の勾配は $g_k=\dfrac m{\Delta t}(2q_k-q_{k-1}-q_{k+1})-\Delta t\,V'(q_k)$。
  $g_k=0$ は $m(q_{k+1}-2q_k+q_{k-1})/\Delta t^2=-V'(q_k)$、すなわち**速度ベルレ法の位置の更新式**
  （第 1 章と接続する）。一様重力では停留点が放物線に**節点上で厳密に**一致する。
- **経路の緩和**: 最急降下法 $q\leftarrow q-\alpha g$、$\alpha=0.9\,\Delta t/(2m)$。
  式にすると拡散方程式の陽解法（拡散係数 $0.45<0.5$ で安定）に等しく、ギザギザが先に消え、全体の形が後から収束する。
  最小点があるときだけ使える（鞍点では発散する）。

### 3.3 ガリレイの相対性（§3, §4）

慣性系 $K'$ が $K$ に対して速度 $\bm V$ で動く: $\bm r=\bm r'+\bm Vt$、$\bm v=\bm v'+\bm V$、$t=t'$。

- 自由粒子: $\mathcal L=\tfrac12mv^2=\tfrac12mv'^2+m\bm V\!\cdot\!\bm v'+\tfrac12mV^2=\mathcal L'+\dfrac{df}{dt}$、
  $f=m\bm V\!\cdot\!\bm r'+\tfrac12mV^2t$。よって $S-S'=f(t_2)-f(t_1)=m\bm V\!\cdot(\bm r'_2-\bm r'_1)+\tfrac12mV^2\tau$。
  **両端を固定した全経路で同じ値**なので、停留経路は変わらない。
- 電車の中の投げ上げ: 電車系 $K'$ で $(0,\ v_0t-\tfrac12gt^2)$、地面系 $K$ で $(Vt,\ v_0t-\tfrac12gt^2)$。
  加速度はどちらも $(0,-g)$（形が同じ）。運動エネルギーは $K$ と $K'$ で異なる（不変ではない）。

### 3.4 質点系（§5）

二体（$m_1,m_2$、万有引力、$G=1$）: $\mathcal L=\tfrac12m_1v_1^2+\tfrac12m_2v_2^2-U(|\bm r_1-\bm r_2|)$、$U=-Gm_1m_2/r$。
EL 方程式は $m_a\ddot{\bm r}_a=-\partial U/\partial\bm r_a$（ニュートンの運動方程式）。
$U$ が相対位置だけに依存するので、全運動量 $\bm P=\sum m_a\bm v_a$ が保存し、重心は等速直線運動をする。
既定の軌道: $m_1=3,\ m_2=1$、近点距離 $r_p=1$、離心率 $e=0.5$ →
$v_{\rm rel}^2=G(m_1+m_2)(1+e)/r_p=6$、長半径 $a=2$、周期 $T=2\pi\sqrt{a^3/(G(m_1+m_2))}=2\pi\sqrt2\approx8.886$。

## 4. 数値手法

| 処理 | 手法 | 刻み・分割 |
|---|---|---|
| 二重振り子 | RK4（`numerics.integrate`） | $h=1$ ms |
| 二体問題 | 速度ベルレ（保存系・ガリレイ変換に共変） | $h=1$ ms |
| 作用の積分 | シンプソン則（`scipy.integrate.simpson`） | 4001 点 |
| 経路の緩和 | 最急降下法（`numerics/optimize.py`） | 節点 41（$\Delta t=0.05$ s）、4000 反復 |

速度ベルレは $K$ と $K'$ で**丸め誤差の範囲で同じ軌道**を出す（加速度が相対位置だけに依存するため、
$\bm V$ を足しても更新式が形を変えない）。これが §3 の不変性の数値的な確認になる。

## 5. パラメータと入出力

| コマンド | 主な引数（既定値） | 出力 |
|---|---|---|
| `least-action` | `--duration 2.0`（重力の例）、`--fps 24` | `outputs/mechanics/least_action_relaxation.mp4`、図 5 枚 |
| `generalized-coordinates` | `--duration 10`、`--theta1-deg 90`、`--theta2-deg 0`、`--fps 30` | `outputs/mechanics/double_pendulum.mp4`、図 2 枚 |
| `galilean` | `--duration 13.3`（$1.5T$）、`--drift-x 0.35`、`--drift-y 0.2`、`--fps 30` | `outputs/mechanics/galilean_two_body.mp4`、図 3 枚 |

図は `outputs/figures/`（各コマンドの `figures()`）。シナリオの数値:

- 重力の例: $m=1$、$g=9.81$、$(t,q)=(0,0)\to(2,0)$。停留経路は頂点 $4.905$ m の放物線、$S_*=-g^2\tau^3/24=-32.08$、$c_1=\tfrac{m\tau}4(\pi/\tau)^2=1.234$。
- 調和振動子の例: $m=1$、$\omega=2\pi$（$\pi/\omega=0.5$ s）、$q_1=1,\ q_2=0$、
  短い $\tau=0.3$（$c_1=+5.26$）と長い $\tau=0.8$（$c_1=-4.81,\ c_2=+4.44$）。

## 6. 検証方法（必須）

| 層 | 検証 | 期待 |
|---|---|---|
| domain（lagrangian） | 3 系の停留経路が EL 方程式（$m\ddot q=-V'$）と端点条件を満たす | 残差 $<10^{-8}$ |
| domain（lagrangian） | 離散作用の勾配が差分と一致する。一様重力で勾配 0 の経路が放物線に厳密に一致 | 相対 $10^{-8}$ / $10^{-12}$ |
| domain（double\_pendulum） | $m_1=m_2,\ l_1=l_2$ の基準振動数が $\sqrt{(2\mp\sqrt2)g/l}$。$m_2\to0$ で単振り子に戻る | 相対 $10^{-12}$ |
| domain（two\_body） | 力 $=-\nabla U$（独立な数値微分） | 相対 $10^{-6}$ |
| application（least\_action） | $S(\varepsilon)$ が 2 次式で、$c_n$ が公式と一致。$\varepsilon=0$ が最小（短い $\tau$）／鞍点（長い $\tau$） | 相対 $10^{-5}$ |
| application（least\_action） | $c_1$ が $\tau=\pi/\omega$ で符号を変える | 零点 $=\pi/\omega$（$10^{-3}$） |
| application（least\_action） | 緩和: 作用が単調に減少し、停留経路（重力は放物線）に収束 | 誤差 $<10^{-6}$ m |
| application（generalized\_coordinates） | 二重振り子の $E$ が保存（RK4）。基準振動の周期が理論値と一致。最大角が $\cos\theta_{\max}$ の式と一致 | $\Delta E/E<10^{-7}$ / $10^{-3}$ |
| application（generalized\_coordinates） | シミュレーションの運動が $\mathcal L$ から数値微分で作った EL 方程式の残差 0 を満たす | 相対 $10^{-4}$ |
| application（galilean） | 作用のずれ $S-S'$ が経路によらず $m V\Delta q'+\tfrac12mV^2\tau$ | 相対 $10^{-8}$ |
| application（galilean） | $K$ で独立に積分した二体の軌道 $=$ $K'$ の軌道 $+\bm Vt$。$\bm P$ 保存、$E$ 保存、重心の等速運動 | $10^{-9}$ / $10^{-12}$ / $10^{-5}$ / $10^{-12}$ |
| 描画 | 煙テスト（`figures()` が全ファイルを書く） | 例外なし |
| 目視 | PDF の該当ページと動画のフレームを `--preview` / `ffmpeg` で確認 | 切れ・重なり・文字化けなし |

## 7. 結果・知見

[knowledge/analytical-mechanics.md](../knowledge/analytical-mechanics.md) に実測値と、ハマった点を書く。

## 8. 未解決・次にやること

- 保存則（§6 エネルギー、§7 運動量、§9 角運動量）とネーターの定理は次の章。
- 経路積分（ファインマン）の見方は読み物として触れるに留め、シミュレーションは未着手。
- 二重振り子のカオス（初期値の敏感性）は未検証。
