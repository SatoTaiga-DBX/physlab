# ドキュメント

このプロジェクトでは「**教科書で理解 → 設計書を書く → 作る → 確かめる → 知見を貯める**」を 1 サイクルとして回す。
コードは `src/`、文書はすべて `docs/` に置く。

| 場所 | 役割 | 書くタイミング |
|---|---|---|
| [textbook/](textbook/README.md) | **教科書**（LaTeX。PDF にして読む）。理論・導出・演習 | 題材の最初（シミュレーションの前後に読む） |
| [design/](design/) | **設計書**。何をどう作るか、どう検証するか | **実装の前**（実装後にずれたら直す） |
| [knowledge/](knowledge/) | **知見**。実装・実行して分かったこと、ハマったこと | 実装・検証の**後** |
| [roadmap.md](roadmap.md) | 取り組む題材の一覧と進捗 | 題材を決めたとき・終えたとき |

教科書は「物理として何が正しいか」、設計書は「それをどう作るか」、知見は「作って何が分かったか」。
同じ内容を繰り返さず、リンクで結ぶ。

## 教科書 (`textbook/`)

| 章 | 題 |
|---|---|
| [第 1 章](textbook/tex/chapters/00-foundations/01-ode-and-numerical-integration.tex) | 微分方程式と数値積分 |
| [第 2 章](textbook/tex/chapters/01-mechanics/01-pendulum.tex) | 単振り子 |
| [第 3 章](textbook/tex/chapters/02-analytical-mechanics/01-equations-of-motion.tex) | 運動方程式（参考: ランダウ＝リフシッツ『力学』§1–§5） |
| [第 4 章](textbook/tex/chapters/02-analytical-mechanics/02-conservation-laws.tex) | 保存法則と力学的相似（参考: 同 §6–§10） |

PDF の作り方・目次・書き方は [textbook/README.md](textbook/README.md)。

## 設計書 (`design/`)

| # | 題名 | 状態 |
|---|---|---|
| [000](design/000-template.md) | テンプレート（コピーして使う） | — |
| [001](design/001-architecture.md) | プロジェクト全体の設計と決定事項 | 実装済み |
| [002](design/002-pendulum.md) | 単振り子 | 実装済み・検証済み |
| [003](design/003-layered-architecture.md) | レイヤー構成（オニオン／クリーンアーキテクチャ流） | 実装済み・検証済み |
| [004](design/004-tex-textbook.md) | 教科書の TeX 化と PDF ビルド（`physlab run textbook`） | 実装済み・検証済み |
| [005](design/005-analytical-mechanics-ch1.md) | 解析力学入門（ランダウ＝リフシッツ第 I 章）のシミュレーション | 実装済み・検証済み |
| [006](design/006-conservation-laws.md) | 保存法則と力学的相似のシミュレーション | 実装済み・検証済み |

## 知見 (`knowledge/`)

| 題名 | 内容 |
|---|---|
| [numerical-integration.md](knowledge/numerical-integration.md) | オイラー / RK4 / ベルレ法の精度とエネルギー保存の違い（実測） |
| [visualization.md](knowledge/visualization.md) | matplotlib で動画を作るときの落とし穴と対処 |
| [architecture.md](knowledge/architecture.md) | 層構造・import-linter・テスト配置を運用して分かったこと |
| [analytical-mechanics.md](knowledge/analytical-mechanics.md) | 最小作用の原理・一般座標・ガリレイ変換のシミュレーションの実測値と、ハマった点 |
| [latex-textbook.md](knowledge/latex-textbook.md) | LuaLaTeX（MiKTeX）で日本語の教科書を組むときの設定・警告・ハマった点 |
| [conservation-laws.md](knowledge/conservation-laws.md) | 対称性と保存量・積分法・力学的相似の実測値と分かったこと |
| [copyright.md](knowledge/copyright.md) | 参考書をもとに教科書を書くときの著作権の調査と方針 |

## 書き方のルール

- 設計書は「検証方法」の節を必須にする（解析解・保存量・極限との照合）。
- 知見には**根拠**（どのシミュレーションで、どのパラメータで、何が分かったか）を添える。
  実測した数値と、推測・未検証のことは区別して書く。
- 教科書の数値は**計算してから**書く。
- 新しい教科書・設計書・知見を足したら、この README の表にも 1 行足す。
