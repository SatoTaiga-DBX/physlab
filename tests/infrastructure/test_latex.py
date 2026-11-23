"""LaTeX ビルド部品: ログの解析と図の存在確認は純粋関数でテストし、実コンパイルは lualatex があるときだけ。"""

import shutil
from pathlib import Path

import pytest

from physlab.infrastructure.latex import (
    LatexBuildError,
    build_pdf,
    error_excerpt,
    find_lualatex,
    find_missing_figures,
    referenced_figures,
    summarize_log,
)

STYLE = Path(__file__).resolve().parents[2] / "docs/textbook/tex/physlab.sty"

LOG_OK = """\
LaTeX Warning: Reference `sec:x' on page 3 undefined on input line 10.
LaTeX Warning: Reference `sec:x' on page 4 undefined on input line 12.
LaTeX Warning: Citation `ll-mech' on page 5 undefined on input line 3.
Overfull \\hbox (67.5pt too wide) in paragraph at lines 7--32
Overfull \\hbox (0.8pt too wide) in paragraph at lines 40--41
LaTeX Warning: Label(s) may have changed. Rerun to get cross-references right.
Output written on book.pdf (20 pages, 443783 bytes).
"""


def test_summarize_log_collects_problems() -> None:
    s = summarize_log(LOG_OK)
    assert s.undefined_references == ("sec:x",)  # 同じ鍵は 1 回だけ
    assert s.undefined_citations == ("ll-mech",)
    assert s.overfull_boxes == (67.5,)  # 1pt に満たないはみ出しは無視する
    assert s.needs_rerun
    assert s.pages == 20


def test_summarize_log_clean() -> None:
    s = summarize_log("Output written on x.pdf (1 page, 100 bytes).")
    assert not s.needs_rerun
    assert s.undefined_references == ()
    assert s.pages == 1


def test_error_excerpt_starts_at_first_error() -> None:
    log = "ok line\n./chapters/a.tex:5: Undefined control sequence.\nl.5 \\foo\nmore\n"
    assert error_excerpt(log).startswith("./chapters/a.tex:5:")


def test_find_missing_figures(tmp_path: Path) -> None:
    tex = tmp_path / "tex"
    tex.mkdir()
    (tex / "a.tex").write_text(
        "\\physfig[0.8]{fig-one}{cap}{lab}\n"
        "\\includegraphics[width=3cm]{fig-two.png}\n"
        "% \\physfig{commented-out}{c}{l}\n",
        encoding="utf-8",
    )
    figures = tmp_path / "figures"
    figures.mkdir()
    (figures / "fig-one.pdf").write_bytes(b"%PDF")
    assert referenced_figures(tex) == {"fig-one", "fig-two.png"}  # コメント行の図は数えない
    assert find_missing_figures(tex, figures) == ["fig-two.png"]


needs_tex = pytest.mark.skipif(find_lualatex() is None, reason="lualatex が無い")


@needs_tex
def test_style_compiles_and_resolves_references(tmp_path: Path) -> None:
    """physlab.sty のすべての環境が通り、問題↔解答・式・章の参照が解決される。"""
    src = tmp_path / "tex"
    src.mkdir()
    shutil.copy(STYLE, src / "physlab.sty")
    (src / "main.tex").write_text(
        r"""\documentclass[a4paper,11pt]{ltjsbook}
\usepackage{physlab}
\begin{document}
\chapter{試験}\label{ch:t}
\chapterinfo{★☆☆}{約 5 分}{2026-01-01}\correspondsto{試験}
\begin{objectives}\begin{enumerate}\item 作用 $S=\int\Lagr\,\dd t$\end{enumerate}\end{objectives}
\begin{defn}{一般座標}{gc}定義 本文\end{defn}
\begin{principle}{原理}{p}原理 本文\end{principle}
\begin{result}{結果}{r}式~\eqref{eq:a} と定義~\ref{def:gc}。\end{result}
\begin{equation}\diff{y}{x}=\pdiff{f}{x}\label{eq:a}\end{equation}
\begin{intuition}x\end{intuition}\begin{note}x\end{note}\begin{simbox}x\end{simbox}
\begin{summary}x\end{summary}
\begin{exercise}[基礎]{a}問題。\end{exercise}
\section{演習の解答}
\begin{solution}{a}答え。第\ref{ch:t}章\end{solution}
\end{document}
""",
        encoding="utf-8",
    )
    result = build_pdf(src, tmp_path / "out", tmp_path / "figs", jobname="t", max_passes=4)
    assert result.pdf.is_file()
    assert result.passes >= 2  # 相互参照の解決に 2 回以上かかる
    assert result.summary.undefined_references == ()
    assert result.summary.pages == 2


@needs_tex
def test_build_failure_reports_log_excerpt(tmp_path: Path) -> None:
    src = tmp_path / "tex"
    src.mkdir()
    (src / "main.tex").write_text(
        "\\documentclass{article}\\begin{document}\\undefinedmacro\\end{document}\n",
        encoding="utf-8",
    )
    with pytest.raises(LatexBuildError, match="undefinedmacro|Undefined control sequence"):
        build_pdf(src, tmp_path / "out", tmp_path / "figs", jobname="bad")
