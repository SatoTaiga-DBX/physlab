"""``physlab`` コマンド: presentation/commands 配下のシミュレーションを探して実行する。

コマンドモジュールの規約は ``physlab.presentation.commands`` を参照。

使い方::

    uv run physlab list
    uv run physlab run pendulum --theta0-deg 150
"""

from __future__ import annotations

import argparse
import importlib
import pkgutil
import sys
from types import ModuleType

from physlab.presentation import commands


def discover() -> dict[str, ModuleType]:
    """commands 配下から規約を満たすモジュールを集め、名前で引ける辞書にする。"""
    found: dict[str, ModuleType] = {}
    for info in pkgutil.walk_packages(commands.__path__, prefix=f"{commands.__name__}."):
        if info.ispkg:
            continue
        module = importlib.import_module(info.name)
        if all(hasattr(module, attr) for attr in ("NAME", "TITLE", "main")):
            found[module.NAME] = module
    return dict(sorted(found.items()))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="physlab", description="物理シミュレーションの実行ランナー"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="実行できるシミュレーションの一覧を表示する")
    run = sub.add_parser("run", help="シミュレーションを実行する")
    run.add_argument("name", help="シミュレーション名（`physlab list` で確認）")
    run.add_argument("args", nargs=argparse.REMAINDER, help="シミュレーション固有の引数")
    ns = parser.parse_args(argv)

    sims = discover()
    if ns.command == "list":
        width = max((len(n) for n in sims), default=0)
        for name, module in sims.items():
            print(f"{name:<{width}}  {module.TITLE}")
        return 0

    if ns.name not in sims:
        print(f"'{ns.name}' というシミュレーションはありません。", file=sys.stderr)
        print(f"候補: {', '.join(sims) or '(なし)'}", file=sys.stderr)
        return 2
    sims[ns.name].main(ns.args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
