# 可視化：matplotlib で動画を作るときの落とし穴

出典: [設計書 002 単振り子](../design/002-pendulum.md)（2026-11-23）。
いずれも数値テストでは見つからず、**動画から 1 フレーム抜いて目で見て**初めて気づいた。

## 動画が実時間より遅く再生される

`stride = round(1 / (fps * dt))` でフレームを間引くと、`1/(30 × 0.001) = 33.3` が 33 に丸まり、
1 コマが 0.033 s（本来は 0.0333 s）進むので約 1% 遅くなる。10 s の指定が 10.13 s の動画になった。

対処: コマ番号 k の時刻 `k/fps` に最も近いステップを使う。

```python
frames = np.round(np.arange(0.0, duration, 1.0 / fps) / dt).astype(int)
```

修正後は 300 フレーム・ちょうど 10.000 s（ffprobe で確認）。

## 表示範囲が動きに足りない

振り子の表示範囲を「支点より下」で固定したため、120° で錘が支点より上に来ると棒が切れた。
範囲は固定値でなく**初期条件から決める**（`y_top = max(0.3, 1.2·(−cos θ0))`）。
パラメータを引数で変えられるようにしたので、固定値はすぐ破綻する。

## 符号付き誤差を symlog に描くと読めない

エネルギー誤差は正負が入れ替わる。symlog 軸では振動する曲線が帯状に潰れて、
ドリフトか有界かが判別しづらかった。**絶対値の対数軸**（`semilogy` 相当）にすると、
「増え続ける / ゆっくりドリフト / 有界に振動」の違いが一目で分かる。
0 に落ちる点があるので `np.maximum(abs(x), 1e-16)` で下駄を履かせる。

## 環境まわり（Windows）

- **日本語フォント**: matplotlib の既定では豆腐（□）になる。`infrastructure/rendering/style.py` で
  Yu Gothic → Meiryo → … の順に、入っているものを使う。`axes.unicode_minus = False` も必要。
- **mp4 出力**: ffmpeg が PATH にあれば `FFMpegWriter`（H.264, `-pix_fmt yuv420p`）。
  `yuv420p` を付けないと一部のプレーヤーで再生できない。幅・高さは偶数にする（1000×500）。
  ffmpeg が無ければ GIF にフォールバックする。
- **Agg バックエンド固定**: ウィンドウを出さずファイルにだけ書くため。
- **日本語と記号の字間**: タイトルに `θ0 = 120°` のように Greek 文字と日本語フォントが混ざると
  字間が崩れることがあった。ラベルは「初期角 120°」のように日本語側に寄せた。

## 確認の手順

```powershell
ffprobe -v error -select_streams v:0 -show_entries stream=width,height,duration,nb_frames -of default=nw=1 outputs\mechanics\pendulum.mp4
ffmpeg -y -loglevel error -ss 1.4 -i outputs\mechanics\pendulum.mp4 -frames:v 1 frame.png
```
