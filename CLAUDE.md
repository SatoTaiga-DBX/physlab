# physlab — プロジェクト規約

物理シミュレーションを Python で作り、動画と教科書で理解するための学習用プロジェクト。
ユーザーは物理の学び直しが目的で、新しい題材を次々に追加していく。
ドキュメントは日本語、コードの識別子は英語。

## コマンド（パッケージ管理は uv のみ。pip は使わない）

```powershell
uv sync                           # 環境構築
uv add <pkg>                      # 依存追加（開発用は uv add --dev <pkg>）
uv run physlab list               # コマンド一覧（シミュレーション・figures・textbook）
uv run physlab run <name> [args]  # 実行（出力は outputs/<topic>/）
uv run physlab run figures        # 教科書に貼る図を作る（outputs/figures/*.pdf）
uv run physlab run textbook       # 教科書（TeX）を PDF にする: 図 → LuaLaTeX → 検査
uv run pytest                     # テスト（層の依存ルールの検査も含む）
uv run lint-imports               # 層の依存ルールだけを検査
uv run ruff check . ; uv run ruff format .
```

PowerShell では `uv run` のビルドログが stderr に出るため、`2>&1` でつなぐと
`NativeCommandError` に見えることがある（失敗ではない。終了コードで判断する）。

## 層構造（オニオン／クリーンアーキテクチャ流）— 詳細は docs/design/003-layered-architecture.md

依存は**外側から内側へだけ**。内側は外側を知らない。

```
presentation ──▶ infrastructure ──▶ application ──▶ domain | numerics
 (CLI 組み立て)    (描画・ファイル)    (ユースケース)    (物理の法則 | 数値解法)
```

| 書くもの | 置き場所 |
|---|---|
| 運動方程式、保存量、解析解、パラメータの dataclass | `src/physlab/domain/<topic>/<name>.py` |
| 積分法・差分スキーム・格子など、物理を知らない計算 | `src/physlab/numerics/` |
| 走らせて結果データ（dataclass）を返す／手法の比較／結果から量を測る | `src/physlab/application/<topic>/<name>.py` |
| matplotlib 描画、動画・画像の保存、出力パス | `src/physlab/infrastructure/` |
| 引数処理、`print`、各層の組み立て（`NAME` / `TITLE` / `main`） | `src/physlab/presentation/commands/<topic>/<name>.py` |
| テスト | `tests/` は src を鏡写し: `tests/<層>/<topic>/test_<name>.py` |

- domain / numerics / application は **matplotlib を import しない**。層違反は `lint-imports`
  （`tests/test_architecture.py`）が検出する。契約は `pyproject.toml` の `[tool.importlinter]`。
- 迷ったら: 「描画やファイルを消しても意味がある計算か？」→ はい: 内側 3 層、いいえ: infrastructure。
  「時間を進める処理か？」→ はい: application、いいえ: domain。
- 描画方法を差し替えたくなるまで、ports（抽象基底クラス）や DI は作らない。

## 開発フロー（必ずこの順）

0. **教科書の章**（TeX。`docs/textbook/tex/chapters/<部>/NN-name.tex`）: 物理の理論・導出・演習。
   `docs/textbook/tex/_template.tex` をコピーして `main.tex` から `\input` する。
   式は途中を省略せず、数値は計算してから書く。問題は `exercise`、解答は章末の `solution`（鍵で相互リンク）。
1. **設計書**: `docs/design/000-template.md` をコピーして `docs/design/NNN-名前.md`。
   支配方程式・仮定・数値手法・**検証方法**を先に決める。
2. **実装**: 内側の層から順に（domain → numerics → application → infrastructure → presentation）。
3. **検証**: 層ごとにテスト。domain は式・極限・特殊値、application は解析解・保存則との照合。
   「動いた」ではなく「物理的に正しい」ことをテストで示す。
4. **目視**: 動画から 1 フレーム抜いて（`ffmpeg -ss <秒> -i x.mp4 -frames:v 1 f.png`）見た目を確認する。
   教科書は `uv run physlab run textbook --preview <ページ>` で PNG にして見る。
   表示範囲の切れ・時間のずれ・文字化け・凡例の重なりは数値テストでは見つからない。
5. **知見の記録**: 分かったこと・ハマったことを `docs/knowledge/` に残す。
   数値は実測値を書き、未検証の推測は「未検証」と明記する。
6. `docs/roadmap.md` のチェックと、`docs/README.md`・`docs/textbook/README.md` の一覧を更新する。
   `uv run physlab run textbook` で未解決の参照が 0 であることを確認する。

設計と実装がずれたら、設計書を直す（コードだけ直して放置しない）。

## 新しい題材を足すチェックリスト

- [ ] 教科書の章 `docs/textbook/tex/chapters/…`（TeX。`main.tex` に `\input` を足す）
- [ ] 設計書 `docs/design/NNN-….md`
- [ ] `domain/<topic>/<name>.py`（＋ `__init__.py` が無ければ追加）
- [ ] `application/<topic>/<name>.py`（結果の dataclass を返す）
- [ ] `infrastructure/rendering/<topic>/<name>.py`（結果を受け取って描く）
- [ ] `presentation/commands/<topic>/<name>.py`（`NAME` / `TITLE` / `main`）→ 自動で `physlab list` に出る
- [ ] 教科書に貼る図があれば、同じコマンドに `figures() -> list[Path]` を足す（`figure_path("名前")` に PDF で出す）。
  `tests/presentation/test_figures.py` の `EXPECTED_FIGURES` にも足す
- [ ] テスト（domain / application は必須、描画は煙テスト、CLI は既存の `test_cli.py` が一覧を検査）
- [ ] 動画のフレームと PDF のページを目視、`uv run pytest` と `uv run ruff check .` を通す
- [ ] 知見・roadmap・各 README の一覧を更新

## 参考書をもとに章を書くとき（著作権。詳細は docs/knowledge/copyright.md）

- このリポジトリは**公開**している。参考書の文章・図・問題を写さない。文を追いかけた言い換えもしない。
- 構成は自分で決め、図と数値は自分のシミュレーションで作り、演習は独自に作る。
- 物理法則・式・導出の筋道は書いてよい。参考にした書名と節番号は出典として示す。
- push の前に上の 3 点と個人情報の混入を確認する。

## 規約

- 単位は SI（m, kg, s, rad）。角度の度数表記は入出力の境界（presentation）でのみ使う。
- 数値積分は `numerics/integrators.py` を使う。手法の比較自体が学習対象なので、
  安易に `solve_ivp` へ置き換えない（必要になったら併用を設計書で決める）。
- matplotlib は Agg 固定（`infrastructure/rendering/style.setup_matplotlib()` を図の作成前に呼ぶ）。
- 動画は実時間と一致させる（k コマ目の時刻 k/fps に最も近い計算ステップを使う）。
  現象が短すぎてスロー再生にする場合は、倍率を**画面に表示**する。
- 教科書の図を作る `figures()` は**引数を取らない**（コマンドの引数を変えても図が別の条件で上書きされないように）。
- テストは `--import-mode=importlib`（同名ファイル可）。テストファイル同士は import できないので、
  共有は `tests/conftest.py` のフィクスチャで行う。テストは `PHYSLAB_OUTPUT_DIR` を tmp に向けるので
  `outputs/` を汚さない。
- `outputs/` は Git 管理外。生成物（動画・図・PDF）をコミットしない。教科書は TeX が原稿で、PDF は作って読む。
- 教科書の TeX は LuaLaTeX + LuaTeX-ja + 原ノ味フォント（MiKTeX）。仕組みは `docs/design/004-tex-textbook.md`。
- matplotlib の数式（mathtext）の文字列は **raw 文字列**（`r"$\theta$"`）にする。普通の文字列だと
  `\t` がタブになり、`$\theta$` が「heta」と出る（実際にやった）。
- テストは `tests/` を src の鏡写しにするが、`tests/docs/`（教科書の TeX ソースの整合性）は層ではない例外。

## Windows / PowerShell の注意

- **シェルのヒアドキュメントにバックスラッシュを含む文章（TeX・正規表現・mathtext）を書かない。**
  ツール経由でバックスラッシュが減ることがある。ファイルは Write / Edit ツールで作る。

- ファイルの編集には Edit / Write ツールを使う。Windows PowerShell 5.1 の
  `Set-Content -Encoding utf8` は **BOM を付ける**、`Get-Content`（`-Encoding` なし）は
  BOM なし UTF-8 を ANSI として読んで**文字化けして見える**。
  どうしても PowerShell で書くときは `[System.IO.File]::WriteAllText(path, text, (New-Object System.Text.UTF8Encoding $false))`。
