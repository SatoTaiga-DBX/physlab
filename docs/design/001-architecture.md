# 001: プロジェクト全体の設計

- 状態: 実装済み（フォルダ構成の部分は [003](003-layered-architecture.md) が置き換えた）
- 作成日: 2026-11-23
- 更新日: 2026-11-23（003 を受けて構成を更新）

## 1. 目的

物理法則をシミュレーションし、動画と教科書で理解する。題材を次々に追加できること、
作るたびに設計と知見が蓄積されることを最優先にする。

## 2. 構成

ソースはオニオン／クリーンアーキテクチャ流の層構造。詳細・理由・置き場所の早見表は
[003 レイヤー構成](003-layered-architecture.md) を参照。

```
src/physlab/
├── domain/          物理の法則（運動方程式・保存量・解析解）
├── numerics/        数値解法（euler / rk4 / velocity_verlet）
├── application/     ユースケース（走らせて結果データを返す）
├── infrastructure/  描画・動画保存・出力先
└── presentation/    `physlab` コマンドと各題材のコマンド
docs/                textbook / design / knowledge / roadmap
tests/               src を鏡写し
```

### シミュレーションの規約

`presentation/commands/<topic>/<name>.py` が `NAME`（実行名）・`TITLE`（一覧表示用）・`main(argv)` を持てば、
`pkgutil.walk_packages` による自動検出で `physlab list` / `physlab run <NAME>` に現れる。
レジストリへの手動登録を無くし、題材追加の手間を最小にするため。

### 共通化の方針

- 共通化するのは**複数の題材で確実に使うもの**だけ（積分器・動画保存・描画設定）。
- 題材固有の物理（運動方程式、エネルギー、解析解）は `domain/<topic>/` に閉じる。
- 偏微分方程式（波動・流体）が始まったら、差分格子・境界条件のヘルパを `numerics/` に足す。
  それまでは作らない（実際の要求が見えてから設計する）。

## 3. 決定事項

| 決定 | 理由 |
|---|---|
| パッケージ管理は uv のみ | ユーザー要望。`pyproject.toml` + `uv.lock` で再現性を確保 |
| Python は 3.13 に固定（`.python-version`） | 3.14 も導入済みだが、科学計算系は新しすぎる版より安定版を選んだ（判断であり、3.14 での不具合を確認したわけではない） |
| `src/` レイアウト + `uv init --package` | `uv run physlab ...` でコマンド化でき、テストが実インストール相当の import で走る |
| 層構造と依存ルール（内向きのみ） | 題材が増えても置き場所に迷わず、物理が描画の都合に引きずられない（[003](003-layered-architecture.md)） |
| 層の依存ルールは import-linter で自動検査 | 約束を守る努力でなく、テストで担保する |
| 積分器は自前実装 | 手法ごとの精度・保存則の違いを比べること自体が学習目的。本格計算には scipy を併用可 |
| 描画は matplotlib（Agg 固定）、動画は ffmpeg の mp4 | 追加の依存が少なく、ffmpeg 導入済みなら高圧縮の mp4 が得られる。無ければ GIF |
| 出力は `outputs/`（Git 管理外） | 動画は再生成できる。環境変数 `PHYSLAB_OUTPUT_DIR` で変更可 |
| lint/format は ruff、テストは pytest（importlib モード） | 標準的で高速。tests を src の鏡写しにするため同名ファイルを許す |
| ドキュメントは日本語、識別子は英語 | 学習ノートとしての読みやすさ |
| 文書は `docs/` に集約（教科書・設計書・知見） | ユーザー要望。コードと文書を分ける |

## 4. 検証方法

- 各題材は解析解または保存量でテストする（`tests/`）。
- 数値解法（`numerics/`）はテスト可能な既知解（調和振動子 y'' = -y）で収束次数を確認する
  （`tests/numerics/test_integrators.py`）。
- 層の依存ルールは `tests/test_architecture.py`（import-linter）で検査する。

## 5. 未解決・次にやること

- 偏微分方程式向けの共通部品（格子・境界条件・可視化）の設計。波動方程式の題材で着手する。
- 対話的に見たい場合の `--show` オプション（現状は Agg 固定でファイル出力のみ）。
- Git 管理（未初期化）。必要になったら `git init`。
