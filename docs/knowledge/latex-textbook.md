# TeX 教科書：LuaLaTeX（MiKTeX）で日本語を組むときの設定とハマった点

出典: [設計書 004](../design/004-tex-textbook.md)（2026-11-23）。
環境: Windows 11、MiKTeX 25.12（winget）、LuaLaTeX + LuaTeX-ja、原ノ味フォント。

## 確認できたこと

- `winget install --id MiKTeX.MiKTeX --scope user --silent` で、**管理者権限なし**で入った。
  実行ファイルは `%LOCALAPPDATA%\Programs\MiKTeX\miktex\bin\x64\`。導入直後に開いていた端末の PATH には載らないので、
  `infrastructure/latex.py` の `find_lualatex()` が既定の場所も探す。
- `initexmf --set-config-value "[MPM]AutoInstall=1"` で、足りないパッケージ（tcolorbox、tikz、siunitx、listings …）が
  初回のビルド時に自動で入った。ビルドのコマンドには MiKTeX のときだけ `--enable-installer` を付ける。
- `\usepackage[haranoaji]{luatexja-preset}` で、原ノ味明朝・ゴシックが MiKTeX から取得され、日本語フォントを
  システムに頼らず組める（他の環境でも同じ見た目になる）。
- ビルド時間: 42 ページ・図 10 枚で **26 秒**（図の生成 11 秒＋LuaLaTeX 2 回。`physlab run textbook`）。
  44 ページ・図 13 枚の現在の版は 2 回のコンパイルで収束する。
- `latexmk.exe` は MiKTeX に入っているが Perl が要るので使わず、ログが「もう一度実行せよ」と言わなくなるまで
  自前で繰り返す（最大 4 回）。
- PDF を画像にして見るには `pdftoppm` が無いので、`pypdfium2`（開発用の依存）を使う。

## ハマったこと

| 症状 | 原因 | 対処 |
|---|---|---|
| `Undefined control sequence: \tcbcolframe` | tcolorbox のスタイルの中で、色を決める前に色を参照した | 色を**引数に取るスタイル**（`physbox=<色名>`）にした |
| 章扉のページ番号が上に出る | `ltjsbook` の `\plainifnotempty` が、fancyhdr のヘッダが空でないため「上にページ番号」の型を選ぶ | `\renewcommand{\plainifnotempty}{\thispagestyle{plain}}` |
| `\headheight is too small (14pt)` | fancyhdr の既定の高さが足りない | geometry で `headheight=18pt` |
| 表がはみ出す（`Overfull \hbox` 67〜141 pt） | 長いパスを `l` 列に入れた | `tabularx` と折り返せる列 `L` |
| 長いコマンドを含む問題文で、かな文字の字間が伸びる | 折り返せない長い文字列があると、行の余りが日本語の字間に配られる | 長いコマンドは `terminal` 環境（コードブロック）に出す |
| `\qty{1e-9}{m}` が「1·10⁻⁹ m」 | `exponent-product={\cdot}` を設定していた | 設定を外して既定の「×」にした |
| `\pdiff[2]{f}{q}` で崩れる | 自作の `\pdiff` は引数 2 つで、省略可能な引数を取らない | `\frac{\partial^2f}{\partial q^2}` と書く |
| `caption` が「Unknown document class」と警告 | `ltjsbook` を知らないだけ | 見た目への影響は出ていないので放置（設計書 004 に記載） |
| `Japanese fonts will be scaled by 0.9247` | 原ノ味と欧文フォントの大きさ合わせ | 無害 |
| matplotlib の PDF で `'created' timestamp seems very low` | fonttools がフォントの日付を読む警告 | 無害（`pdf.fonttype=42` で TrueType を埋め込むときに出る） |

## 設計上のうまくいった点

- **問題と解答を鍵（`key`）で結ぶ**ことで、「解答は p.N →」「問 N の解答」が相互リンクになり、
  対応が崩れると**未解決の参照でビルドが失敗**する。`tests/docs/test_textbook_sources.py` が TeX を走らせずに先に検出する。
- **図を `figures()` で作って `outputs/figures/` に PDF で出す**ことで、数値を直すと図も再生成され、
  本文の数値との食い違いに気づける。図の参照がなければコンパイル前に失敗する。
- 本文に書く数値は、実行して**出力された値**を使った。演習の解答に載せる数値（$\tau=1$ の $S_*$、$c$、ドリフトを変えたときの差）は、
  書く前に実行して確かめ、予想と違う箇所（ドリフト $(0,0)$ の差は「$10^{-16}$ 程度」でなく**厳密に 0**）を直した。

## 未検証

- TeX Live・Linux・macOS での動作。
- 日本語フォントを変えたときの `scaled` と字面の違い。
- タグ付き PDF（アクセシビリティ）。
