"""PDF のページ指定の解釈と、PNG への書き出し。"""

from pathlib import Path

from matplotlib import pyplot as plt

from physlab.infrastructure.pdf_preview import count_pages, parse_pages, render_pages


def test_parse_pages_handles_lists_ranges_and_all() -> None:
    assert parse_pages("1,3,5-7", 10) == [1, 3, 5, 6, 7]
    assert parse_pages("all", 3) == [1, 2, 3]
    assert parse_pages("2,2,9-12", 10) == [2, 9, 10]  # 重複と範囲外は捨てる


def test_render_pages_writes_png(tmp_path: Path) -> None:
    pdf = tmp_path / "sample.pdf"
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    fig.savefig(pdf)
    plt.close(fig)

    assert count_pages(pdf) == 1
    (png,) = render_pages(pdf, [1], tmp_path / "png", scale=0.5)
    assert png.name == "sample-p001.png"
    assert png.stat().st_size > 0
