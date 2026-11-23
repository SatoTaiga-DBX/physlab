# 第1-1章 単振り子

- 難易度: ★★☆
- 所要時間: 約 50 分
- 更新日: 2026-11-23

## 学習目標

1. 単振り子の運動方程式 $\ddot\theta = -(g/L)\sin\theta$ を導ける。
2. 小角近似で単振動になること、周期が振幅によらない（等時性）ことを説明できる。
3. エネルギー保存から、厳密な周期が楕円積分で書けることを導ける。
4. 振幅が大きいと周期が伸びる理由と、その大きさを見積もれる。
5. 位相空間で単振り子の運動を読める。

## 前提

- 微積分、三角関数、ニュートンの運動方程式。
- 数値積分（§8）は [第0-1章](../00-foundations/01-ode-and-numerical-integration.md) を先に読むと分かりやすい。

## 1. 導入

ガリレオはピサの大聖堂のランプの揺れから、「振れ幅が変わっても周期は同じ」と気づいたと伝えられる。
これは**近似として**正しい。振れ幅が小さければ周期はほとんど変わらないが、大きく振らせると、はっきり遅くなる。

この章の問い: **振れ幅が大きくなると、周期はどれだけ、なぜ伸びるのか？**

## 2. 設定

質量 $m$ の質点を、長さ $L$ の軽くて伸びない棒で吊るす。支点まわりに摩擦はなく、重力加速度は $g$。
鉛直下向きからの角度を $\theta$（反時計回りを正）とする。質点の位置は

$$
(x, y) = (L\sin\theta,\ -L\cos\theta).
$$

自由度は $\theta$ の 1 つ。

## 3. 運動方程式

質点にはたらく力は重力 $m g$（鉛直下向き）と棒の張力 $S$（棒に沿って支点向き）。
重力を、棒の**接線方向**（$\theta$ が増える向き）と**棒に沿った方向**に分ける。

- 接線方向の成分: $-mg\sin\theta$（$\theta$ を減らす向き＝復元力）
- 棒に沿った成分: $mg\cos\theta$（張力と釣り合う側）

質点は半径 $L$ の円周上を動くので、接線方向の加速度は $L\ddot\theta$。接線方向のニュートンの運動方程式は

$$
mL\ddot\theta = -mg\sin\theta
\quad\Longrightarrow\quad
\boxed{\ddot\theta = -\frac{g}{L}\sin\theta}.
$$

質量 $m$ が消えた。**周期は質量によらない**。

棒に沿った方向は、向心加速度 $L\dot\theta^2$ を与える：$S - mg\cos\theta = mL\dot\theta^2$。
つまり張力は $S = m(g\cos\theta + L\dot\theta^2)$。運動方程式には効かない（棒が伸びないという拘束で決まる）。

## 4. 小角近似と単振動

$|\theta| \ll 1$ では $\sin\theta = \theta - \theta^3/6 + \cdots \approx \theta$。方程式は

$$
\ddot\theta = -\omega_0^2\,\theta, \qquad \omega_0 = \sqrt{g/L}
$$

となり、一般解は $\theta(t) = A\cos\omega_0 t + B\sin\omega_0 t$。静かに放した場合（$\theta(0)=\theta_0$、$\dot\theta(0)=0$）は

$$
\theta(t) = \theta_0\cos\omega_0 t, \qquad T_0 = \frac{2\pi}{\omega_0} = 2\pi\sqrt{\frac{L}{g}}.
$$

周期 $T_0$ は振幅 $\theta_0$ を含まない。これが**等時性**。$L = 1\,\mathrm{m}$ なら $T_0 = 2.0061\,\mathrm{s}$。

### 近似がどれだけ粗いか

小角近似は、復元力 $\sin\theta$ を $\theta$ で置き換えている。$\theta_0 = 120° = 2.094\,\mathrm{rad}$ では
$\sin\theta_0 = 0.866$ だが $\theta_0 = 2.094$ で、**近似は復元力を 2.4 倍に過大評価**する。
本当の振り子は近似より弱い力で引き戻されるので、**ゆっくり**動く＝周期が伸びる。
動画で、橙（近似）が青（厳密）より先に一周するのはこのため。

## 5. エネルギー保存

運動方程式に $\dot\theta$ を掛けると

$$
\dot\theta\,\ddot\theta = -\frac gL\sin\theta\;\dot\theta
\quad\Longleftrightarrow\quad
\frac{d}{dt}\left[\tfrac12\dot\theta^2 - \frac gL\cos\theta\right] = 0.
$$

$mL^2$ を掛けて整えると、保存量（力学的エネルギー、最下点を基準）が得られる：

$$
E = \tfrac12 mL^2\dot\theta^2 + mgL(1-\cos\theta) = \text{一定}.
$$

静かに放した場合 $E = mgL(1-\cos\theta_0)$ なので、角速度は

$$
\dot\theta^2 = \frac{2g}{L}\left(\cos\theta - \cos\theta_0\right). \tag{5.1}
$$

最下点（$\theta=0$）の速さは $v = L\dot\theta = \sqrt{2gL(1-\cos\theta_0)}$。

## 6. 厳密な周期

(5.1) を $\dot\theta = d\theta/dt$ とみて変数分離する。振り子は $\theta = 0 \to \theta_0$ に $T/4$ かけて動くから

$$
\frac T4 = \int_0^{\theta_0}\frac{d\theta}{\dot\theta}
= \sqrt{\frac{L}{2g}}\int_0^{\theta_0}\frac{d\theta}{\sqrt{\cos\theta-\cos\theta_0}}.
$$

倍角の公式 $\cos\theta = 1 - 2\sin^2(\theta/2)$ で書き換えると
$\cos\theta-\cos\theta_0 = 2\left[\sin^2(\theta_0/2) - \sin^2(\theta/2)\right]$。
$k = \sin(\theta_0/2)$ とおき、置換 $\sin(\theta/2) = k\sin\varphi$（$\varphi: 0\to\pi/2$）を行う。

- 根号の中: $\cos\theta-\cos\theta_0 = 2k^2\cos^2\varphi$
- 微分: $\tfrac12\cos(\theta/2)\,d\theta = k\cos\varphi\,d\varphi$、かつ $\cos(\theta/2) = \sqrt{1-k^2\sin^2\varphi}$
  なので $d\theta = \dfrac{2k\cos\varphi}{\sqrt{1-k^2\sin^2\varphi}}\,d\varphi$

よって被積分関数は $\dfrac{d\theta}{\sqrt{\cos\theta-\cos\theta_0}} = \dfrac{\sqrt2\,d\varphi}{\sqrt{1-k^2\sin^2\varphi}}$。$\cos\varphi$ が約分されて、

$$
\boxed{T = 4\sqrt{\frac Lg}\;K(k) = \frac{4}{\omega_0}K\!\left(\sin\frac{\theta_0}{2}\right)},
\qquad
K(k) = \int_0^{\pi/2}\frac{d\varphi}{\sqrt{1-k^2\sin^2\varphi}}.
$$

$K(k)$ は**第 1 種完全楕円積分**。$k\to0$ で $K\to\pi/2$ となり、$T\to 2\pi/\omega_0 = T_0$（小角近似）に戻る。
（SciPy の `ellipk(m)` は $m = k^2$ を引数に取る。コードの `period_exact` はこの式そのもの。）

## 7. 周期は振幅でどう変わるか

$K(k)$ を $k$ で展開すると $K = \tfrac\pi2\left[1 + \tfrac14k^2 + \tfrac9{64}k^4 + \cdots\right]$。
$k^2 = \sin^2(\theta_0/2) \approx \theta_0^2/4$ を使うと

$$
\frac{T}{T_0} = 1 + \frac{\theta_0^2}{16} + \frac{11\,\theta_0^4}{3072} + \cdots
$$

（$\theta_0$ は rad）。振幅が小さいうちは**周期の伸びは振幅の 2 乗に比例**する。

計算値（$T/T_0$ の厳密値と、第 1 項まで取った近似 $1+\theta_0^2/16$）:

| 初期角 $\theta_0$ | $T/T_0$（厳密） | $1+\theta_0^2/16$ | 近似の誤差 |
|---|---|---|---|
| 10° | 1.0019 | 1.0019 | 0.00% |
| 30° | 1.0174 | 1.0171 | −0.03% |
| 60° | 1.0732 | 1.0685 | −0.43% |
| 90° | 1.1803 | 1.1542 | −2.21% |
| 120° | 1.3729 | 1.2742 | −7.19% |
| 150° | 1.7622 | 1.4284 | −18.94% |
| 170° | 2.4394 | 1.5502 | −36.45% |
| 175° | 2.8777 | 1.5831 | −44.99% |
| 179° | 3.9011 | 1.6100 | −58.73% |

- 30° まで振っても周期の伸びは 1.7% しかない。ガリレオの観察は実用上正しい。
- 120° で 37% 伸びる。第 1 項だけの近似は 7% ずれる。
- $\theta_0\to180°$（真上）で周期は**発散**する。$k\to1$ で $K(k)\simeq\ln\dfrac{4}{\sqrt{1-k^2}}$ なので

$$
T \simeq \frac{4}{\omega_0}\ln\frac{4}{\cos(\theta_0/2)} \qquad(\theta_0\to\pi).
$$

  179° で式の値 $T/T_0 = \dfrac{2}{\pi}\ln\dfrac{4}{\cos 89.5°} = 3.901$ は表の厳密値と一致する。
  物理的には、真上は不安定な釣り合い点で、そこへ向かう振り子は頂点付近でほとんど止まり、通過に長い時間がかかる。

## 8. 位相空間で見る

位置 $\theta$ と速度 $\omega=\dot\theta$ を軸にした平面（位相空間）では、エネルギー保存 (5.1) より運動は
**等エネルギー曲線**の上を動く：

$$
\omega = \pm\frac1L\sqrt{2\left[E/m - gL(1-\cos\theta)\right]}.
$$

- 小さい振幅: ほぼ**楕円**（単振動）。
- 振幅が大きくなる: 楕円が $\theta$ 方向に引き伸ばされ、「魚の目」の形になる。
- $E/m = 2gL$（真上で静止）: **分離線**（セパラトリクス）。これより小さい $E$ は振動、大きい $E$ は回転運動。
  本章は振動（$\theta_0<\pi$）だけを扱う。

閉じた曲線＝周期運動。数値解法がうまいかどうかは、**この閉曲線から外れないか**で見られる。
前進オイラー法は外へ螺旋状に広がり（エネルギーが増える）、速度ベルレ法は閉曲線の近くに留まる。
これが [第0-1章](../00-foundations/01-ode-and-numerical-integration.md) の話。

## 9. シミュレーションで確かめる

```powershell
uv run physlab run pendulum                      # 初期角 120°（既定）
uv run physlab run pendulum --theta0-deg 10      # 小角: 近似とほぼ重なる
uv run physlab run pendulum --theta0-deg 170     # 大振幅: 周期が約 2.4 倍
```

出力は `outputs/mechanics/pendulum.mp4`。見るところ:

- 左: 厳密解（青）と小角近似（橙）の 2 つの振り子。初めは重なっていて、だんだん位相がずれる。
- 右: $\theta(t)$。橙のほうが山の間隔が短い。タイトルに 2 つの周期が出る。

初期角 120°、$L = 1\,\mathrm{m}$ での値:

| 量 | 理論値 | 備考 |
|---|---|---|
| 小角近似の周期 $T_0$ | 2.0061 s | $2\pi\sqrt{L/g}$ |
| 厳密な周期 $T$ | 2.7541 s | $\tfrac{4}{\omega_0}K(\sin 60°)$、$T/T_0=1.3729$ |

数値解（RK4、$h=1\,\mathrm{ms}$）から求めた周期は、$\theta_0 = 10, 60, 120, 170°$ のすべてで厳密な周期と
相対誤差 $10^{-4}$ 以内で一致する（`tests/application/mechanics/test_pendulum.py`）。
つまり「楕円積分の公式」と「微分方程式を数値で解いた結果」は、独立に同じ答えに到達する。

## 演習

**問1（基礎）** 長さ $L=0.25\,\mathrm{m}$ の振り子の、小角での周期 $T_0$ を求めよ。

<details>
<summary>解答</summary>

$T_0 = 2\pi\sqrt{0.25/9.81} = 2\pi\times0.1596 \approx 1.003\,\mathrm{s}$。
（長さを 1/4 にすると周期は 1/2。$T_0\propto\sqrt L$。）

</details>

**問2（基礎）** $\theta_0 = 120°$、$L=1\,\mathrm{m}$ で放した振り子の、最下点での速さを求めよ。

<details>
<summary>解答</summary>

§5 より $v=\sqrt{2gL(1-\cos\theta_0)} = \sqrt{2\times9.81\times1\times1.5} = \sqrt{29.43} \approx 5.42\,\mathrm{m/s}$。
$\cos120° = -0.5$ なので $1-\cos\theta_0 = 1.5$。質量は不要。

</details>

**問3（応用）** $\theta_0=60°$ で、近似式 $T\approx T_0(1+\theta_0^2/16)$ の値と厳密値を比べよ。
2 項目 $\tfrac{11}{3072}\theta_0^4$ まで入れるとどうなるか。

<details>
<summary>解答</summary>

$\theta_0 = 1.0472\,\mathrm{rad}$。第 1 項まで: $1 + 1.0966/16 = 1.0685$。厳密値は 1.0732（誤差 −0.43%）。
$11\theta_0^4/3072 = 11\times1.2026/3072 = 0.00431$ を足すと $1.0728$。厳密値 1.0732 との差は 0.04% まで縮まる。

</details>

**問4（応用）** $\theta_0\to\pi$ で周期が発散する理由を、(a) 物理的に (b) 式の上で説明せよ。

<details>
<summary>解答</summary>

(a) 真上は不安定な釣り合い点（力が 0）。そこに近い高さから放すと、頂点付近で速度がほとんど 0 になり、
通過に非常に長い時間がかかる。
(b) $k=\sin(\theta_0/2)\to1$ で $K(k)=\int_0^{\pi/2}d\varphi/\sqrt{1-k^2\sin^2\varphi}$ の被積分関数が
$\varphi=\pi/2$ で発散する（積分値は対数的に発散）。§7 の $T\simeq(4/\omega_0)\ln[4/\cos(\theta_0/2)]$。

</details>

**問5（応用）** 初期角速度が 0 でないとき（$\dot\theta(0)=\omega_i$）、小角近似の解は？ また、真上を通り抜ける
（回転運動になる）ために必要な、最下点での最小の速さは？

<details>
<summary>解答</summary>

小角近似: $\theta(t)=\theta_0\cos\omega_0t+(\omega_i/\omega_0)\sin\omega_0t$（コードの `small_angle_solution`）。
回転になる条件: $E/m>2gL$。最下点で $E/m=\tfrac12v^2$ なので $v>2\sqrt{gL}$。$L=1\,\mathrm{m}$ で $v>6.26\,\mathrm{m/s}$。

</details>

**問6（シミュレーション）** `--theta0-deg 10` と `--theta0-deg 170` の動画を見比べ、
(a) 橙と青のずれ方の違い、(b) 右のグラフで 10 秒間に何周するか、を答えよ。表の値と合っているか確かめよ。

<details>
<summary>解答の目安</summary>

10° ではほぼ重なったまま（$T/T_0=1.0019$）。170° では $T = 2.4394\,T_0 = 4.89\,\mathrm{s}$ なので、10 秒間に約 2 周。
橙（$T_0 = 2.006\,\mathrm{s}$）は約 5 周する。

</details>

## まとめ

- 運動方程式は $\ddot\theta=-(g/L)\sin\theta$。質量に依存しない。
- 小角近似では単振動で周期 $T_0=2\pi\sqrt{L/g}$。振幅によらない（等時性）。
- 厳密な周期は楕円積分 $T=(4/\omega_0)K(\sin(\theta_0/2))$。振幅が大きいほど伸びる。
- 周期の伸びは小振幅では $\theta_0^2/16$。120° で 37%、真上に近づくと発散。
- 位相空間では等エネルギー曲線の上を動く。数値解法の良し悪しはその曲線から外れないかで見える。

## 発展

- **ラグランジュ形式**: $\mathcal L=\tfrac12mL^2\dot\theta^2+mgL\cos\theta$ に
  オイラー＝ラグランジュ方程式を使うと、張力を考えずに同じ運動方程式が得られる。
- **二重振り子**: 自由度が 2 になり、カオスが現れる。
- **減衰・強制振動**: 摩擦と外力を加え、共鳴を調べる。
- **回転運動**: 分離線を越えて回り続ける領域。

## 対応するファイル

| 種類 | 場所 |
|---|---|
| 設計書 | [002 単振り子](../../design/002-pendulum.md) |
| 物理（domain） | `src/physlab/domain/mechanics/pendulum.py`（`deriv`, `energy`, `period_exact`, `small_angle_solution`） |
| ユースケース | `src/physlab/application/mechanics/pendulum.py`（`simulate`, `compare_with_small_angle`） |
| 描画 | `src/physlab/infrastructure/rendering/mechanics/pendulum.py` |
| テスト | `tests/domain/mechanics/test_pendulum.py`, `tests/application/mechanics/test_pendulum.py` |
| 知見 | [numerical-integration.md](../../knowledge/numerical-integration.md), [visualization.md](../../knowledge/visualization.md) |

## 参考文献

- J. R. Taylor『Classical Mechanics』（振動・振り子の章）。
- ランダウ＝リフシッツ『力学』（楕円積分による周期の導出）。
- ファインマン『ファインマン物理学 I』（調和振動、振り子）。
