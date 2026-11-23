# ドキュメント

このプロジェクトでは「**教科書で理解 → 設計書を書く → 作る → 確かめる → 知見を貯める**」を 1 サイクルとして回す。
コードは `src/`、文書はすべて `docs/` に置く。

| 場所 | 役割 | 書くタイミング |
|---|---|---|
| [textbook/](textbook/README.md) | **教科書**。理論・導出・演習。読んで学ぶ | 題材の最初（シミュレーションの前後に読む） |
| [design/](design/) | **設計書**。何をどう作るか、どう検証するか | **実装の前**（実装後にずれたら直す） |
| [knowledge/](knowledge/) | **知見**。実装・実行して分かったこと、ハマったこと | 実装・検証の**後** |
| [roadmap.md](roadmap.md) | 取り組む題材の一覧と進捗 | 題材を決めたとき・終えたとき |

教科書は「物理として何が正しいか」、設計書は「それをどう作るか」、知見は「作って何が分かったか」。
同じ内容を繰り返さず、リンクで結ぶ。

## 教科書 (`textbook/`)

| 章 | 題 |
|---|---|
| [0-1](textbook/00-foundations/01-ode-and-numerical-integration.md) | 微分方程式と数値積分 |
| [1-1](textbook/01-mechanics/01-pendulum.md) | 単振り子 |

目次と書き方は [textbook/README.md](textbook/README.md)。

## 設計書 (`design/`)

| # | 題名 | 状態 |
|---|---|---|
| [000](design/000-template.md) | テンプレート（コピーして使う） | — |
| [001](design/001-architecture.md) | プロジェクト全体の設計と決定事項 | 実装済み |
| [002](design/002-pendulum.md) | 単振り子 | 実装済み・検証済み |
| [003](design/003-layered-architecture.md) | レイヤー構成（オニオン／クリーンアーキテクチャ流） | 実装済み・検証済み |

## 知見 (`knowledge/`)

| 題名 | 内容 |
|---|---|
| [numerical-integration.md](knowledge/numerical-integration.md) | オイラー / RK4 / ベルレ法の精度とエネルギー保存の違い（実測） |
| [visualization.md](knowledge/visualization.md) | matplotlib で動画を作るときの落とし穴と対処 |
| [architecture.md](knowledge/architecture.md) | 層構造・import-linter・テスト配置を運用して分かったこと |

## 書き方のルール

- 設計書は「検証方法」の節を必須にする（解析解・保存量・極限との照合）。
- 知見には**根拠**（どのシミュレーションで、どのパラメータで、何が分かったか）を添える。
  実測した数値と、推測・未検証のことは区別して書く。
- 教科書の数値は**計算してから**書く。
- 新しい教科書・設計書・知見を足したら、この README の表にも 1 行足す。
