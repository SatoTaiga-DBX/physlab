"""PDF の指定ページを PNG にする（目視確認用）。"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path


def parse_pages(spec: str, n_pages: int) -> list[int]:
    """``"1,3,5-7"`` や ``"all"`` を、1 始まりのページ番号の列にする（範囲外は捨てる）。"""
    pages: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if part == "all":
            pages += range(1, n_pages + 1)
        elif "-" in part:
            start, end = part.split("-", 1)
            pages += range(int(start), int(end) + 1)
        elif part:
            pages.append(int(part))
    return [p for p in dict.fromkeys(pages) if 1 <= p <= n_pages]


def render_pages(pdf: Path, pages: Iterable[int], out_dir: Path, scale: float = 1.4) -> list[Path]:
    """``pages``（1 始まり）を ``out_dir/<pdf名>-pNNN.png`` に書き出し、パスの一覧を返す。"""
    import pypdfium2 as pdfium  # 遅延 import: 目視確認のときだけ必要（開発用の依存）

    out_dir.mkdir(parents=True, exist_ok=True)
    document = pdfium.PdfDocument(str(pdf))
    written: list[Path] = []
    for p in pages:
        image = document[p - 1].render(scale=scale).to_pil()
        path = out_dir / f"{pdf.stem}-p{p:03d}.png"
        image.save(path)
        written.append(path)
    return written


def count_pages(pdf: Path) -> int:
    import pypdfium2 as pdfium

    return len(pdfium.PdfDocument(str(pdf)))
