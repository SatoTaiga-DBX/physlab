# アーキテクチャ：層構造・import-linter・テスト配置を運用して分かったこと

出典: [設計書 003 レイヤー構成](../design/003-layered-architecture.md)（実施日 2026-11-23）

## 実測・確認したこと

1. **移行は結果を変えなかった。** 単振り子（120°）を移行前後で実行し、周期（2.7541 / 2.0061 s）、
   エネルギー誤差（オイラー +3.21e+00、RK4 −7.05e-08、ベルレ最大 1.84e-04）、動画（300 フレーム・10.000 s）が
   完全に一致。動画の 1 フレームも見た目が同じ。
2. **import-linter は意図的な違反を検出した。** 初回は 34 ファイル・42 依存で 2 契約とも成立。
   わざと壊して確認:
   - `application` から `infrastructure` を import → `layers` 契約が BROKEN。違反した import が 1 行で示される。
   - `domain` に `import matplotlib` → `forbidden` 契約が BROKEN。直接の違反に加え、
     `application → domain → matplotlib` という**間接経路**も報告される。
3. **pytest の `--import-mode=importlib` で、同名の `test_pendulum.py` を 3 つの層に置けた。**
   全 30 件が成功（domain 8 / numerics 5 / application 9 / infrastructure 3 / presentation 4 / architecture 1）。
4. **import-linter は pytest から呼べる。** `importlinter.cli.lint_imports(config_filename=...)` が
   契約成立で 0 を返すので、`tests/test_architecture.py` に入れて `uv run pytest` だけで層の検査も走る。

## 設計上のうまくいった点

- **application が結果を dataclass で返し、描画は受け取るだけ**にしたため、描画のテストが
  「小さな結果を作って渡し、ファイルができるか見る」だけで済んだ。計算を再実行しなくてよい。
- **描画関数が出力パスを引数で受け取る**（パスの決定は presentation）ため、テストは `tmp_path` を渡せる。
  描画が出力先の規約を知らなくてよい。
- **domain のテストに時間積分が要らない。** 式・極限・特殊値（最下点・水平・真上の加速度、120° で
  T/T0 = 1.3729 など）だけで書け、速くて失敗原因が分かりやすい。
- `tests/conftest.py` の autouse フィクスチャで出力先を tmp に向けたので、テストが `outputs/` を汚さない。

## コスト（トレードオフ）

- 1 題材 = 約 4 つのソースファイル（domain / application / rendering / command）＋ 層ごとのテスト。
  型にはまる代わりに、小さな題材でもファイルが増える。チェックリスト（CLAUDE.md）で手順を固定して対処。
- ports（抽象基底クラス）や DI は作っていない。描画の差し替え要求が出るまで不要と判断（未検証の判断）。

## ハマったこと（Windows / PowerShell）

- **PowerShell 5.1 の `Set-Content -Encoding utf8` は BOM を付ける。** 検証のために一時的に
  ソースを書き換えて復元した際、2 ファイルの先頭に BOM が残った。
  `Get-Content`（`-Encoding` なし）は BOM なし UTF-8 を ANSI として読み、表示が文字化けして見える
  （ファイルは壊れていない）。ファイル編集は Edit/Write ツールを使うか、
  `[System.IO.File]::WriteAllText(path, text, (New-Object System.Text.UTF8Encoding $false))` を使う。
- 検証のつもりで入れた「違反」がシェルのエスケープ（`\n`）の誤りで構文エラーになり、
  **検査が通ったように見えない**まま終わりかけた。「壊してみて検出されること」の確認では、
  壊した結果が**想定した種類の違反**として報告されたかまで見る。

## 未検証・今後

- 偏微分方程式（波動・流体）でも、格子・ステンシルを `numerics/` に置く構成が自然かは未検証。波動の題材で確かめる。
- 題材が増えたら、トピック間（`mechanics` ↔ `waves`）の import を import-linter の `independence` 契約で禁止するか検討。
- 層構造のコストが題材の規模に見合うかは、題材を数個足してから再評価する。
