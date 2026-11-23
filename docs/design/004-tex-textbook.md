# 004: 教科書の TeX 化と PDF ビルド

- 状態: 検証済み
- 作成日: 2026-11-23
- 実装:
  - スタイル・本体: [docs/textbook/tex/](../textbook/tex/)（`physlab.sty`、`main.tex`、`chapters/`）
  - ビルド部品: [infrastructure/latex.py](../../src/physlab/infrastructure/latex.py)、
    [infrastructure/pdf_preview.py](../../src/physlab/infrastructure/pdf_preview.py)
  - コマンド: [presentation/commands/textbook.py](../../src/physlab/presentation/commands/textbook.py)、
    [presentation/commands/figures.py](../../src/physlab/presentation/commands/figures.py)
- 実行: `uv run physlab run textbook`

## 1. 目的

Markdown の教科書は、式は出るが、**問題から解答への参照・章をまたぐ参照・表・図の配置**が読みづらかった
（折りたたみ `<details>` はプレビューによって開閉できない）。
TeX にして PDF で読めるようにし、次の 3 つを実現する。

- 式・表・図・参照が**書籍の体裁**で出る。問題と解答は**相互にリンク**する（ページ番号つき）。
- 図は Python のシミュレーションから**ベクタ画像**で貼り込む（手でスクリーンショットを撮らない）。
- ビルドは**1 コマンド**（図の生成 → PDF → 検査）。

## 2. 方針と選定理由

| 決めたこと | 理由 |
|---|---|
| 原稿の正は **TeX**。Markdown の章は廃止 | 2 重管理（MD と TeX）は必ずずれる。MD の章は git の履歴に残る |
| エンジンは **LuaLaTeX + LuaTeX-ja** | 日本語組版が現行の標準的な選択肢。OpenType フォントを直接使える |
| フォントは **原ノ味**（`luatexja-preset` の `haranoaji`） | TeX の配布物に含まれ、環境依存しない（システムのフォントを使うと他の人が組めない） |
| 処理系は **MiKTeX**（winget） | 必要なパッケージだけを初回に自動取得でき、導入が軽い（TeX Live でも動くはずだが**未検証**） |
| 図は **matplotlib → PDF**（データ由来）と **TikZ**（概念図） | 数値から作る図は再現可能に、模式図は TeX の中で式と同じ書体にする |
| ビルドは Python のコマンドから `lualatex` を呼ぶ | `latexmk` は Perl が要る。自前の反復（最大 4 回）は短く、ログも検査できる |

## 3. ディレクトリ

```
docs/textbook/
  README.md                  読み方・ビルド方法・章の一覧
  tex/
    main.tex                 本体（章を並べる）
    physlab.sty              共通スタイル（環境・マクロ）
    chapters/
      00-foundations/NN-名前.tex     部ごとにフォルダを切る
      01-mechanics/…  02-analytical-mechanics/…
      appendix/symbols.tex           記号表
      bibliography.tex               参考文献
    _template.tex            新しい章の雛形
outputs/figures/*.pdf        図（生成物。Git 管理外）
outputs/textbook/            PDF とログ（生成物。Git 管理外）
```

## 4. 章の書き方（`physlab.sty` の環境）

| 用途 | 書き方 |
|---|---|
| 章の頭 | `\chapter{題}\label{ch:…}` → `\chapterinfo{難易度}{所要時間}{更新日}` → `\correspondsto{…}` → `objectives` 環境 |
| 定義・原理・結果 | `\begin{defn}{題}{鍵}` / `principle` / `result`（章ごとに番号。`\ref{def:鍵}` `prin:` `res:`） |
| 補足の囲み | `intuition`（イメージ）、`note`（注）、`simbox`（シミュレーションで確かめる）、`summary` |
| 問題と解答 | 本文に `\begin{exercise}[種別]{鍵}`、章末の `\section{演習の解答}` に `\begin{solution}{鍵}` |
| 図 | `\physfig[幅]{ファイル名}{キャプション}{ラベル}`（`outputs/figures/` から読む） |
| コード | コマンドは `terminal` 環境、識別子は `\code{}`、パスは `\path{}`（スラッシュで折り返せる） |
| 単位 | `\qty{9.81}{m/s^2}`、`\ang{120}`（siunitx） |

**問題と解答の対応**: 鍵（`key`）で結ぶ。問題の末尾に「解答は p.N →」、解答の見出しに「問 N の解答」が出て、
どちらもクリックで行き来できる。鍵が対応しないとビルドが「未解決の参照」で失敗し、
`tests/docs/test_textbook_sources.py` が先に検出する。

## 5. 図のパイプライン

1. 各シミュレーションのコマンドに任意の `figures() -> list[Path]` を持たせる
   （[commands/\_\_init\_\_.py](../../src/physlab/presentation/commands/__init__.py) の規約）。
   動画は含めず、数秒で終わる静止画だけを `outputs/figures/<名前>.pdf` に出す。
2. `physlab run figures` が `figures()` を持つ全コマンドを呼ぶ（名前を指定すれば一部だけ）。
3. `physlab run textbook` は、図の生成 → 参照先の図の存在確認 → LuaLaTeX（`\physlabfigdir` に図の場所を渡す）→ 検査。

## 6. ビルドと検査

`build_pdf()` は、ログが「もう一度実行せよ」と言わなくなるまで `lualatex` を繰り返す（最大 4 回）。
ログから次を拾い、`physlab run textbook` が報告する。

| 検査 | 結果 |
|---|---|
| TeX のエラー | 失敗。ログの該当箇所を表示 |
| 未解決の `\ref` / `\cite` | 失敗（PDF は作られるが終了コード 1） |
| 10 pt 以上はみ出した行・表 | 警告（件数と最大値） |
| 参照している図が無い | コンパイル前に失敗 |

`--preview 3,5-7` で、指定ページを PNG にして `outputs/textbook/preview/` に出す（目視確認用）。

## 7. 検証

| 何を | どうやって |
|---|---|
| ログの解析（未解決参照・はみ出し・再実行） | `tests/infrastructure/test_latex.py`（TeX 不要の純粋関数） |
| スタイルの全環境が通り、参照が解決される | 同ファイルの実コンパイル試験（`lualatex` があるときだけ） |
| 問題と解答の対応・ファイルの存在・ラベルの重複 | `tests/docs/test_textbook_sources.py` |
| PDF → PNG | `tests/infrastructure/test_pdf_preview.py` |
| 見た目 | `--preview` で PNG にして目視 |

## 8. リポジトリに入れるもの・入れないもの

| 入れる | 入れない（`.gitignore`） |
|---|---|
| `src/`, `tests/`, `docs/`（TeX の原稿・設計書・知見） | `outputs/`（動画・図・PDF）、`.venv/`、各種キャッシュ |
| `pyproject.toml`, `uv.lock`, `.python-version` | TeX の中間ファイル（`*.aux` `*.log` `*.toc` …） |

PDF は生成物なので入れない。読むときは `uv run physlab run textbook` で作る。

## 9. 未検証・既知の制約

- TeX Live・Linux・macOS では**未検証**（MiKTeX on Windows のみ）。
- 日本語フォントは原ノ味に固定。変えるときは `physlab.sty` の `luatexja-preset` のオプションを変える。
- `caption` パッケージが `ltjsbook` を「未知のクラス」と警告する（見た目に影響は出ていない）。
