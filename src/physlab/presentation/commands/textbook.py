"""`physlab run textbook`: 教科書（TeX）を PDF にする。図の生成 → LuaLaTeX → 検査の順。

設計書: docs/design/004-tex-textbook.md
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from physlab.infrastructure.latex import (
    LatexBuildError,
    LatexNotFoundError,
    build_pdf,
    find_missing_figures,
)
from physlab.infrastructure.pdf_preview import count_pages, parse_pages, render_pages
from physlab.infrastructure.storage import FIGURES_TOPIC, output_dir
from physlab.presentation.commands.figures import generate_all

NAME = "textbook"
TITLE = "教科書（docs/textbook/tex）を PDF にビルドする（図の生成を含む）"
TOPIC = "textbook"
DEFAULT_TEX_DIR = "docs/textbook/tex"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument(
        "--tex-dir", default=DEFAULT_TEX_DIR, help=f"TeX ソースの場所（既定 {DEFAULT_TEX_DIR}）"
    )
    p.add_argument("--no-figures", action="store_true", help="図を作り直さず、前回のものを使う")
    p.add_argument(
        "--preview", metavar="PAGES", help="PDF の指定ページを PNG にする（例: 3,5-7 / all）"
    )
    p.add_argument("--scale", type=float, default=1.4, help="--preview の拡大率")
    args = p.parse_args(argv)

    tex_dir = Path(args.tex_dir)
    if not (tex_dir / "main.tex").is_file():
        print(f"{tex_dir / 'main.tex'} がありません。", file=sys.stderr)
        return 2

    if not args.no_figures:
        generate_all()
    figures_dir = output_dir(FIGURES_TOPIC)
    missing = find_missing_figures(tex_dir, figures_dir)
    if missing:
        print(
            f"TeX が参照する図が {figures_dir} にありません: {', '.join(missing)}", file=sys.stderr
        )
        print("図を作る関数 figures() を持つコマンドを確認してください。", file=sys.stderr)
        return 1

    print("LuaLaTeX でビルド中...")
    try:
        result = build_pdf(tex_dir, output_dir(TOPIC), figures_dir)
    except (LatexNotFoundError, LatexBuildError) as error:
        print(error, file=sys.stderr)
        return 1

    s = result.summary
    print(f"  -> {result.pdf}（{s.pages} ページ、{result.passes} 回のコンパイル）")
    if s.overfull_boxes:
        count, worst = len(s.overfull_boxes), max(s.overfull_boxes)
        print(f"  注意: はみ出した行・表が {count} か所（最大 {worst:.0f}pt）")

    if args.preview:
        pages = parse_pages(args.preview, count_pages(result.pdf))
        for png in render_pages(result.pdf, pages, output_dir(TOPIC) / "preview", args.scale):
            print(f"  -> {png}")

    broken = [*s.undefined_references, *s.undefined_citations]
    if broken:
        print(
            f"未解決の参照があります（\\ref / \\cite の鍵）: {', '.join(broken)}", file=sys.stderr
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
