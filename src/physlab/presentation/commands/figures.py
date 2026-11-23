"""`physlab run figures`: 教科書に貼る図を、``figures()`` を持つ全コマンドから作る。

設計書: docs/design/004-tex-textbook.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

from physlab.presentation.registry import figure_providers

NAME = "figures"
TITLE = "教科書に貼る図を作る（outputs/figures/。名前を指定すると、そのコマンド分だけ）"


def generate_all(only: list[str] | None = None) -> list[Path]:
    """``figures()`` を持つコマンドを順に呼び、作った図のパスを返す。"""
    providers = figure_providers()
    unknown = [n for n in (only or []) if n not in providers]
    if unknown:
        raise ValueError(f"図を作るコマンドではありません: {', '.join(unknown)}（候補: {', '.join(providers)}）")
    made: list[Path] = []
    for name, module in providers.items():
        if only and name not in only:
            continue
        print(f"[{name}] 図を作成中...")
        paths = module.figures()
        for path in paths:
            print(f"  -> {path}")
        made += paths
    return made


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("names", nargs="*", help="図を作るコマンド名（省略すると全部）")
    args = p.parse_args(argv)
    try:
        made = generate_all(args.names or None)
    except ValueError as error:
        print(error)
        return 2
    print(f"{len(made)} 個の図を作りました。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
