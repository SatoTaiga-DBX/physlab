# physlab

物理法則を Python でシミュレーションし、**動画で眺めて・教科書で理解する**ための学習用ラボ。
力学の基本から始めて、波動・流体へと発展させていく。

- 言語 / 管理: Python 3.13 + [uv](https://docs.astral.sh/uv/)
- 計算: NumPy / SciPy、描画・動画: matplotlib（ffmpeg があれば mp4、なければ GIF）
- 進め方: **教科書 → 設計書 → 実装 → 検証 → 知見**（詳細は [docs/README.md](docs/README.md)）

## クイックスタート

```powershell
uv sync                          # 依存関係をインストール
uv run physlab list              # 実行できるコマンド一覧
uv run physlab run pendulum      # 単振り子を実行 → outputs/mechanics/ に動画と図
uv run physlab run pendulum --theta0-deg 30   # 引数でパラメータを変えられる
uv run physlab run least-action  # 最小作用の原理（作用の地形、経路が真の運動に収束する動画）
uv run physlab run textbook      # 教科書を PDF にする → outputs/textbook/physlab-textbook.pdf
uv run pytest                    # テスト（層の依存ルールの検査を含む）
```

出力（動画・図・PDF）は `outputs/` に保存される。`outputs/` は Git 管理外（再生成できるため）。
教科書は **LaTeX** で書いてあり、PDF にして読む（初回だけ `winget install MiKTeX.MiKTeX`）。
読み方と目次は [教科書](docs/textbook/README.md)。最初は第 1 章・第 2 章、解析力学は第 3 章から
（ランダウ＝リフシッツ『力学』第 I 章に沿って、最小作用の原理などを図とシミュレーションで見る）。

## ディレクトリ構成

```
physics/
├── docs/                         文書（コードとは分ける）
│   ├── textbook/                 教科書: 理論・導出・演習（tex/ が原稿。PDF は作って読む）
│   ├── design/                   設計書: 何をどう作るか（実装の前に書く）
│   ├── knowledge/                知見: やって分かったこと
│   └── roadmap.md                取り組む題材の一覧
├── src/physlab/                  ソース（オニオン／クリーンアーキテクチャ流の層構造）
│   ├── domain/                   ② 物理の法則（運動方程式・保存量・解析解）  純粋・描画なし
│   ├── numerics/                 ② 数値解法（積分器など）                    純粋・物理を知らない
│   ├── application/              ③ ユースケース（走らせて結果データを返す）
│   ├── infrastructure/           ④ 外界との接点（matplotlib 描画・ffmpeg・出力先）
│   └── presentation/             ⑤ 入口（`physlab` コマンドと各題材のコマンド）
├── tests/                        src を鏡写し（tests/<層>/<トピック>/test_<題材>.py）
└── outputs/                      生成物（Git 管理外）
```

依存は**外側から内側へだけ**:
`presentation → infrastructure → application → domain | numerics`。
内側の層は matplotlib を知らない。守られているかは `uv run lint-imports`（と `uv run pytest`）が検査する。
各層の下には題材のトピック（`mechanics/`, 将来 `waves/`, `fluids/`）が同じ名前で並ぶ。
設計の理由は [設計書 003](docs/design/003-layered-architecture.md)。

## 新しい題材を追加する

チェックリストは [CLAUDE.md](CLAUDE.md#新しい題材を足すチェックリスト)。要点:

1. 教科書の章と設計書を書く（`docs/textbook/tex/_template.tex`, `docs/design/000-template.md`）
2. 内側の層から実装する: `domain` → `application` → `infrastructure` → `presentation`
   （`presentation/commands/<トピック>/<名前>.py` に `NAME` / `TITLE` / `main(argv)` があれば自動で一覧に出る）
3. 層ごとにテストを書く（解析解・保存則との照合）
4. 動かして動画を確認し、分かったことを `docs/knowledge/` に記録する
