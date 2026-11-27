# 教科書

物理を学び直すための**読み物**。式の導出・直観・数値での確かめ方・演習を 1 章ずつにまとめる。
シミュレーションで「見る」前後に、ここで「理解する」ためのもの。

原稿は **LaTeX（LuaLaTeX + LuaTeX-ja）**で、PDF にして読む。Markdown 版は廃止した
（式・参照・表・図を書籍の体裁で出すため。理由は [設計書 004](../design/004-tex-textbook.md)）。

## PDF を作って読む

```powershell
winget install MiKTeX.MiKTeX       # 初回だけ。TeX 処理系（足りないパッケージは初回のビルド時に自動で入る）
uv sync
uv run physlab run textbook        # 図の生成 → PDF → 検査。outputs/textbook/physlab-textbook.pdf ができる
```

- 図だけ作り直す: `uv run physlab run figures`
- 図は作り直さず TeX だけ組み直す: `uv run physlab run textbook --no-figures`
- 指定ページを画像にして確認する: `uv run physlab run textbook --no-figures --preview 3,5-7`
  （`outputs/textbook/preview/` に PNG が出る）
- 問題と解答は相互にリンクしている。問題の末尾の「解答は p.N →」、解答の見出し「問 N の解答」をクリックする。

## 目次

| 章 | 題 | 原稿 | 対応するシミュレーション | 状態 |
|---|---|---|---|---|
| **第 I 部 準備** | | | | |
| 第 1 章 | 微分方程式と数値積分 | [tex](tex/chapters/00-foundations/01-ode-and-numerical-integration.tex) | `pendulum`（積分法の比較） | 執筆済み |
| **第 II 部 力学の入口** | | | | |
| 第 2 章 | 単振り子 | [tex](tex/chapters/01-mechanics/01-pendulum.tex) | `pendulum` | 執筆済み |
| **第 III 部 解析力学**（参考: ランダウ＝リフシッツ『力学』） | | | | |
| 第 3 章 | 運動方程式（参考: §1 一般座標 〜 §5 質点系のラグランジアン） | [tex](tex/chapters/02-analytical-mechanics/01-equations-of-motion.tex) | `generalized-coordinates`、`least-action`、`galilean` | 執筆済み |
| 第 4 章 | 保存法則（参考: §6 エネルギー 〜 §10 力学的相似） | [tex](tex/chapters/02-analytical-mechanics/02-conservation-laws.tex) | `symmetry`、`center-of-mass`、`areal-velocity`、`similarity` | 執筆済み |
| 第 5 章 | 衝突と散乱（参考: §16 粒子の崩壊 〜 §20 微小角散乱） | [tex](tex/chapters/02-analytical-mechanics/03-collisions.tex) | `decay`、`elastic-collision`、`scattering`、`rutherford`、`small-angle` | 執筆済み |
| （未定） | 運動方程式の積分（1 次元運動、中心力場、ケプラー問題） | | （予定） | 未着手 |
| （未定） | 微小振動 | | （予定） | 未着手 |
| （未定） | 正準方程式 | | （予定） | 未着手 |
| **波動・流体** | （予定: 弦の振動、干渉、回折、拡散、移流、ナビエ–ストークス） | | | 未着手 |

章番号は `main.tex` の並び順で自動的に付く。「運動方程式の積分」は参考書では衝突の前（第 III 章）なので、
書いたら衝突の章の前に入れる予定（そのとき衝突の章は第 6 章になる）。衝突の章は、必要な中心力の積分を章の中で導いている。
題材の候補は [../roadmap.md](../roadmap.md)。

> **著作物について**: 解析力学の部は、ランダウ＝リフシッツ『力学』が扱う話題を参考にした**学習ノート**で、
> 同書の解説書ではない。同書の文章・図・問題を写して作ったものではなく、図は本プロジェクトのコードで作り、
> 演習は独自に用意している。原稿は生成 AI（Claude）の支援で書いており、原典との照合はしていない。
> 考え方と根拠は [../knowledge/copyright.md](../knowledge/copyright.md)、書誌は `tex/chapters/bibliography.tex`。

## 読み方

1. 各章の「この章で学ぶこと」を読み、前提の章を確認する。
2. 本文は**式を自分で導出しながら**読む（途中式を省略しない書き方にしてある）。
3. 「シミュレーションで確かめる」のコマンドを実際に動かし、動画と数値を見る。
4. 演習を解いてから、章末の解答（リンクで飛べる）を見る。

## 章の書き方（新しい章を足すとき）

[tex/\_template.tex](tex/_template.tex) をコピーして `tex/chapters/<部>/NN-名前.tex` に置き、
[tex/main.tex](tex/main.tex) の `\part` の下に `\input` する。環境・マクロは
[tex/physlab.sty](tex/physlab.sty) と [設計書 004](../design/004-tex-textbook.md) の「章の書き方」を見る。

| 節 | 内容 |
|---|---|
| この章で学ぶこと / 前提 | `objectives` 環境。先に読む章は `\ref{ch:…}` で参照 |
| 導入 | 何が不思議か、なぜ知りたいか（直観から入る） |
| 本論 | 設定 → 方程式 → 解く → 結果の意味。**導出は飛ばさない** |
| シミュレーションで確かめる | `simbox` 環境にコマンド、**実測した数値**と理論値の比較、図は `\physfig` |
| 演習 | `exercise` 環境（基礎・応用・シミュレーション）。解答は章末の `solution` 環境（鍵で対応づける） |
| まとめ / 発展 | `summary` 環境で要点 3〜5 個。次に学ぶこと |

ルール:

- 数値は**実際に計算して**から書く（暗算や記憶で書かない）。理論値と実測値は区別する。
- 記号は付録の記号表（`tex/chapters/appendix/symbols.tex`）にそろえる。足したら表にも足す。
- 設計書（何をどう作るか）や知見（実装して分かったこと）の内容は繰り返さず、`\path{…}` でファイルを示す。
- TeX を書いたら `uv run physlab run textbook` を通し、未解決の参照が 0 であることを確認する。
