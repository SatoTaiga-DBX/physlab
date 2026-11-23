"""``physlab`` コマンド: presentation/commands 配下のシミュレーションを探して実行する。

コマンドモジュールの規約は ``physlab.presentation.commands`` を参照。

使い方::

    uv run physlab list
    uv run physlab run pendulum --theta0-deg 150
    uv run physlab run textbook          # 教科書の PDF をビルドする
"""

from __future__ import annotations

import argparse
import sys

from physlab.presentation.registry import discover

__all__ = ["discover", "main"]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="physlab", description="物理シミュレーションの実行ランナー"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="実行できるコマンドの一覧を表示する")
    run = sub.add_parser("run", help="コマンド（シミュレーション・教科書ビルドなど）を実行する")
    run.add_argument("name", help="コマンド名（`physlab list` で確認）")
    run.add_argument("args", nargs=argparse.REMAINDER, help="コマンド固有の引数")
    ns = parser.parse_args(argv)

    sims = discover()
    if ns.command == "list":
        width = max((len(n) for n in sims), default=0)
        for name, module in sims.items():
            print(f"{name:<{width}}  {module.TITLE}")
        return 0

    if ns.name not in sims:
        print(f"'{ns.name}' というコマンドはありません。", file=sys.stderr)
        print(f"候補: {', '.join(sims) or '(なし)'}", file=sys.stderr)
        return 2
    result = sims[ns.name].main(ns.args)
    return result if isinstance(result, int) else 0  # main が整数を返したら終了コードにする


if __name__ == "__main__":
    raise SystemExit(main())
