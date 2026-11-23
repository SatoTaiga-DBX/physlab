# 003: レイヤー構成（オニオン／クリーンアーキテクチャ流のフォルダ設計）

- 状態: 実装済み
- 作成日: 2026-11-23
- 関連: [001 全体設計](001-architecture.md)（本書がフォルダ構成の部分を置き換える）

## 1. 目的

題材（力学 → 波動 → 流体…）が増えてソースが膨らんでも、

- 「どのコードをどこに置くか」で迷わない
- 物理の法則（変えたくない核）が、描画やファイル出力の都合に引きずられない
- テストの置き場所が一意に決まる
- 構成の約束が守られているかを機械的に確かめられる

ようにする。

## 2. 調査結果（要点）

- クリーン / オニオン / ヘキサゴナルは、**ソースコードの依存を内側（ビジネスロジック側）にだけ向ける**点で共通。
  違いは強調点で、ヘキサゴナルは「内と外」の境界 1 つ、オニオンは内側にも複数の輪を作る**層**重視、
  クリーンは**ユースケース**を第一級にして依存ルールを明文化したもの。
  ([Milan Jovanović の比較](https://milanjovanovic.tech/blog/clean-architecture-vs-onion-vs-hexagonal),
  [Thoughtworks](https://www.thoughtworks.com/insights/blog/architecture/demystify-software-architecture-patterns))
- Python の実例は `domain / application / infrastructure / presentation` の 4 層が多いが、
  **標準の配置は無く**、小規模では重くしすぎないことが推奨される
  ([clean-python](https://github.com/nens/clean-python) ほか)。
- 科学計算ソフトでも「数値エンジン」「物理モデル」「可視化」を分離するのが定石
  (例: [soliton_solver](https://arxiv.org/pdf/2603.24370) は solver / theories / visualization を分けている)。
- 依存ルールは [import-linter](https://import-linter.readthedocs.io/en/stable/get_started/configure/) の
  `layers`（上位から下位へのみ許可）と `forbidden`（外部パッケージの禁止も可）契約で自動検査できる。
- pytest は `--import-mode=importlib` にするとテストのファイル名が重複してよい
  （代償としてテストファイル同士の import ができなくなる）。
  ([pytest ドキュメント](https://docs.pytest.org/en/latest/pythonpath.html))

## 3. 方針：何を採り、何を採らないか

| 採る | 採らない（今は） |
|---|---|
| 依存は内向きのみ（外側が内側を import する） | ports の抽象基底クラス / Protocol の先回り定義 |
| 層ごとにフォルダを分ける | DI コンテナ、Repository パターン |
| 物理法則を最内層に置き、描画・ファイル・CLI を外へ追い出す | 1 題材ごとに全層へクラスを作る儀式 |
| import-linter で依存ルールを自動検査 | |

採らない理由: DB や外部 API が無く、実装を差し替える要求も今は無いので、抽象化は「使われない間接層」になる。
**描画方法を差し替えたくなった**（例: 対話表示、別ライブラリ）時点で、application 側に `Protocol` を置いて逆転させる。

## 4. 層の定義

オニオンの輪との対応と、依存してよい先。上の層ほど外側。

| 層 | フォルダ | 責務 | import してよい先 |
|---|---|---|---|
| ⑤ presentation（最外） | `presentation/` | CLI。引数を解釈し、application と infrastructure をつなぐ（組み立て役） | すべて |
| ④ infrastructure | `infrastructure/` | 外界との接点: matplotlib 描画、ffmpeg、ファイル出力先 | application, domain, numerics |
| ③ application | `application/` | ユースケース。物理＋数値解法を組み合わせ、**結果データ**（dataclass）を返す | domain, numerics |
| ② domain（最内） | `domain/` | 物理の法則: パラメータ、運動方程式、エネルギー、解析解。**純粋な関数とデータのみ** | （なし） |
| ② numerics（最内） | `numerics/` | 数値解法: 積分器など。物理を知らない汎用の純計算 | （なし） |

- オニオンの Domain Services に当たる `numerics/` は、物理に依存しないので `domain/` と**並列の最内層**にした
  （互いに import しない）。
- 内側の 3 層（domain / numerics / application）は **matplotlib を import しない**（描画を知らない）。

```
presentation ──▶ infrastructure ──▶ application ──▶ domain
      │                                  │            numerics
      └──────────────▶ （どの内側にも直接依存してよい）
```

## 5. フォルダ構成

各層の下に、題材のトピック（`mechanics/`, `waves/`, `fluids/`）を同じ名前で切る。

```
src/physlab/
├── domain/
│   └── mechanics/pendulum.py          PendulumParams / 運動方程式 / energy / period_exact / 小角近似解
├── numerics/
│   └── integrators.py                 euler / rk4 / velocity_verlet
├── application/
│   └── mechanics/pendulum.py          simulate / compare_with_small_angle / compare_energy_drift（結果 dataclass）
├── infrastructure/
│   ├── storage.py                     output_path（出力先の規約）
│   └── rendering/
│       ├── style.py                   matplotlib 設定（Agg・日本語フォント）
│       ├── writers.py                 save_animation / save_figure
│       └── mechanics/pendulum.py      動画・図の描画（結果データを受け取って描くだけ）
└── presentation/
    ├── cli.py                         `physlab` コマンド本体（コマンドの自動検出）
    └── commands/mechanics/pendulum.py NAME / TITLE / main(argv)

tests/                                 src と同じ階層を鏡写しにする
├── conftest.py                        出力先を tmp に差し替える（outputs/ を汚さない）
├── domain/mechanics/test_pendulum.py
├── numerics/test_integrators.py
├── application/mechanics/test_pendulum.py
├── infrastructure/rendering/test_pendulum.py   （煙テスト）
├── presentation/test_cli.py
└── test_architecture.py               import-linter の契約を pytest から検査
```

## 6. 置き場所の早見表（「この関数はどこ？」）

| 書くもの | 置き場所 |
|---|---|
| 運動方程式、保存量（エネルギー等）、解析解、パラメータの dataclass | `domain/<topic>/<name>.py` |
| 題材をまたいで使う物理の補助（定数、単位換算など） | `domain/shared/`（必要になったら作る） |
| 積分法・差分スキーム・格子・境界条件など物理を知らない計算 | `numerics/` |
| 「このパラメータで走らせて結果を返す」「手法を比較する」「結果から周期を測る」 | `application/<topic>/<name>.py` |
| グラフ・動画の描画、フォント設定、保存、出力パス | `infrastructure/` |
| コマンドライン引数、`print`、各層の組み立て | `presentation/commands/<topic>/<name>.py` |

迷ったら: **「描画やファイルを消しても意味がある計算か？」** → はい: domain / numerics / application、いいえ: infrastructure。
さらに **「時間を進める（シミュレーションを走らせる）処理か？」** → はい: application、いいえ: domain。

## 7. テストの置き方

- `tests/` は `src/physlab/` と同じ階層を**鏡写し**にする。壊れた層がすぐ分かる。
- 種類を層で分ける:
  - domain: 法則そのものの確認（極限・厳密解の公式・エネルギー式）。数値積分を使わない。
  - numerics: 既知解での収束次数・保存則。
  - application: 時間積分した結果を解析解・保存則と照合（周期、エネルギー保存）。
  - infrastructure: 煙テスト（ファイルができる）。見た目は目視で確認（CLAUDE.md の開発フロー 4）。
  - presentation: CLI の入出力。
- `pyproject.toml` で `--import-mode=importlib` を指定し、同名ファイル（各層の `test_pendulum.py`）を許す。
  テストファイル同士は import できないので、共有は `conftest.py` のフィクスチャで行う。

## 8. 機械的な検査（import-linter）

`pyproject.toml` に 2 つの契約を置き、`uv run lint-imports` と `tests/test_architecture.py` で検査する。

1. **layers**: `presentation > infrastructure > application > (domain | numerics)`。下位から上位への import を禁止。
   `domain` と `numerics` は互いも import 禁止。
2. **forbidden**: domain / numerics / application は `matplotlib` を import しない。

## 9. トレードオフと対策

- **1 題材 = 約 5 ファイル**（各層に 1 つ）になる。「どの層に何を書くか」が型にはまる利点と引き換え。
  → 追加手順は [CLAUDE.md](../../CLAUDE.md) / [README](../../README.md) のチェックリストに書いた。
- 題材間（mechanics と waves）の import は今は禁止していない。肥大化したら
  import-linter の `independence` 契約でトピック間を独立させる。

## 10. 移行（旧 → 新）

| 旧 | 新 |
|---|---|
| `core/integrators.py` | `numerics/integrators.py` |
| `core/style.py` | `infrastructure/rendering/style.py` |
| `core/animation.py` の `save_animation`, `save_figure` | `infrastructure/rendering/writers.py` |
| `core/animation.py` の `output_path` | `infrastructure/storage.py` |
| `simulations/mechanics/pendulum.py`（物理） | `domain/mechanics/pendulum.py` |
| 同（`simulate`, 周期測定, エネルギー比較） | `application/mechanics/pendulum.py` |
| 同（`make_video`, `make_energy_figure`） | `infrastructure/rendering/mechanics/pendulum.py` |
| 同（`main`, 引数処理） | `presentation/commands/mechanics/pendulum.py` |
| `cli.py` | `presentation/cli.py` |

## 11. 検証方法

- 移行前後で**同じ数値**になること（回帰）。移行前の値: 周期 厳密 2.7541 s / 小角近似 2.0061 s、
  エネルギー最終誤差 オイラー +3.21e+00・RK4 −7.05e-08、最大絶対誤差 ベルレ 1.84e-04、動画 300 フレーム 10.000 s。
- 既存テスト（12 件）が層別に再配置されても全て通ること。
- `uv run lint-imports` が通ること。**意図的に違反を入れて検出されること**も確認する。

## 12. 結果・知見

[knowledge/architecture.md](../knowledge/architecture.md) に記録。

## 13. 未解決・次にやること

- 偏微分方程式（波動・流体）で格子・ステンシル・境界条件を `numerics/` にどう置くか（波動の題材で設計）。
- 描画を差し替えたくなったら、application 側に出力用の `Protocol` を置く（今は作らない）。
