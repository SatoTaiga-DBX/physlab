# NNN: <題材名>

- 状態: 下書き / 実装中 / 実装済み / 検証済み
- 作成日: YYYY-MM-DD
- 教科書: `docs/textbook/tex/chapters/<部>/NN-<name>.tex`
- 実装（層ごと）:
  - `src/physlab/domain/<topic>/<name>.py`
  - `src/physlab/application/<topic>/<name>.py`
  - `src/physlab/infrastructure/rendering/<topic>/<name>.py`
  - `src/physlab/presentation/commands/<topic>/<name>.py`
- 実行: `uv run physlab run <name>`

## 1. 目的

何を理解したいか。この動画・図で「何が見えれば成功か」を 1〜3 行で。

## 2. 物理モデル

- 支配方程式（できれば式の導出の要点も）
- 仮定・近似（無視しているもの：摩擦、粘性、…）
- 変数と単位（SI）

## 3. 数値手法

- 離散化・積分法とその選定理由
- 刻み幅 / 格子幅、安定性条件（CFL 条件など）
- 想定される誤差の種類

## 4. パラメータと入出力

| 名前 | 既定値 | 意味 |
|---|---|---|
| | | |

出力: `outputs/<topic>/<name>.mp4` など。

## 5. 可視化の設計

動画のレイアウト、何を同時に見せるか、色・軸の決め方。
教科書に貼る静止画は、コマンドの `figures()` で `outputs/figures/` に PDF で出す（名前も決めておく）。

## 6. 検証方法（必須）

「物理的に正しい」ことをどう確かめるか。実装前に決める。

- 解析解との比較:
- 保存量（エネルギー・運動量・質量など）:
- 極限（小さい/大きいパラメータ）での既知の振る舞い:
- 対応するテスト: 層ごとに `tests/domain/<topic>/test_<name>.py`（式・極限・特殊値）、
  `tests/application/<topic>/test_<name>.py`（解析解・保存則との照合）

## 7. 結果・知見

実装後に埋める。分かったことは `docs/knowledge/` に切り出してリンクする。

## 8. 未解決・次にやること
