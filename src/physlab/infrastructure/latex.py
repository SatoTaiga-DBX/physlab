"""LuaLaTeX で教科書の PDF を作る。外部コマンド ``lualatex`` を呼ぶだけで、物理は知らない。

設計書: docs/design/004-tex-textbook.md
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


class LatexNotFoundError(RuntimeError):
    """lualatex が見つからない。"""


class LatexBuildError(RuntimeError):
    """TeX のコンパイルが失敗した（メッセージにログの該当箇所を含む）。"""


@dataclass(frozen=True)
class LogSummary:
    """ログから読み取った、ビルド後に人が見るべき点。"""

    needs_rerun: bool
    undefined_references: tuple[str, ...]
    undefined_citations: tuple[str, ...]
    overfull_boxes: tuple[float, ...]  # はみ出した幅 [pt]（大きいものだけ）
    pages: int | None


@dataclass(frozen=True)
class BuildResult:
    pdf: Path
    log: Path
    passes: int
    summary: LogSummary


# ---- lualatex の場所 -------------------------------------------------------------------


def find_lualatex() -> str | None:
    """PATH から lualatex を探す。無ければ MiKTeX の既定のインストール先も見る。

    MiKTeX をユーザーインストールした直後は、開いている端末の PATH にまだ載っていない。
    """
    found = shutil.which("lualatex")
    if found:
        return found
    candidates: list[Path] = []
    if local := os.environ.get("LOCALAPPDATA"):
        candidates.append(Path(local) / "Programs/MiKTeX/miktex/bin/x64/lualatex.exe")
    if program_files := os.environ.get("ProgramFiles"):
        candidates.append(Path(program_files) / "MiKTeX/miktex/bin/x64/lualatex.exe")
    return next((str(c) for c in candidates if c.exists()), None)


# ---- ログの解析（TeX を実行せずにテストできる純粋関数） ----------------------------------

_UNDEF_REF = re.compile(r"LaTeX Warning: Reference `([^']+)' on page \d+ undefined")
_UNDEF_CITE = re.compile(r"LaTeX Warning: Citation `([^']+)' on page \d+ undefined")
_OVERFULL = re.compile(r"Overfull \\hbox \(([\d.]+)pt too wide\)")
_PAGES = re.compile(r"Output written on .*? \((\d+) pages?")
_RERUN = re.compile(r"Rerun to get|may have changed|has changed\. Rerun")
_OVERFULL_THRESHOLD_PT = 10.0


def summarize_log(text: str) -> LogSummary:
    """lualatex のログから、未解決の参照・はみ出し・再実行の要否を拾う。"""
    pages = _PAGES.search(text)
    return LogSummary(
        needs_rerun=bool(_RERUN.search(text)),
        undefined_references=tuple(dict.fromkeys(_UNDEF_REF.findall(text))),
        undefined_citations=tuple(dict.fromkeys(_UNDEF_CITE.findall(text))),
        overfull_boxes=tuple(
            w for w in (float(x) for x in _OVERFULL.findall(text)) if w >= _OVERFULL_THRESHOLD_PT
        ),
        pages=int(pages.group(1)) if pages else None,
    )


def error_excerpt(text: str, context: int = 12) -> str:
    """ログの最初の ``!`` エラー行とその周辺を取り出す（無ければ末尾）。"""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("!") or re.match(r"^\S+\.tex:\d+:", line):
            return "\n".join(lines[i : i + context])
    return "\n".join(lines[-context:])


# ---- 図の存在確認 ----------------------------------------------------------------------

_FIGURE_USE = re.compile(r"\\(?:physfig|includegraphics)(?:\[[^\]]*\])?\{([^}]+)\}")
_FIGURE_SUFFIXES = ("", ".pdf", ".png", ".jpg", ".jpeg")


def referenced_figures(tex_dir: Path) -> set[str]:
    """TeX ソースが ``\\physfig`` / ``\\includegraphics`` で参照している図の名前。"""
    names: set[str] = set()
    for tex in tex_dir.rglob("*.tex"):
        if tex.name.startswith("_"):  # 雛形（_template.tex）は対象外
            continue
        text = re.sub(r"(?m)^\s*%.*$", "", tex.read_text(encoding="utf-8"))  # コメント行は除く
        names.update(_FIGURE_USE.findall(text))
    return names


def find_missing_figures(tex_dir: Path, figures_dir: Path) -> list[str]:
    """参照されているのに ``figures_dir`` に無い図の名前（ソート済み）。"""
    return sorted(
        name
        for name in referenced_figures(tex_dir)
        if not any((figures_dir / f"{name}{s}").is_file() for s in _FIGURE_SUFFIXES)
    )


# ---- ビルド ----------------------------------------------------------------------------


def build_pdf(
    tex_dir: Path,
    out_dir: Path,
    figures_dir: Path,
    *,
    main: str = "main.tex",
    jobname: str = "physlab-textbook",
    max_passes: int = 4,
) -> BuildResult:
    """``tex_dir/main`` をコンパイルして ``out_dir/<jobname>.pdf`` を作る。

    目次・相互参照のため、ログが「もう一度実行せよ」と言わなくなるまで繰り返す
    （最大 ``max_passes`` 回）。
    """
    exe = find_lualatex()
    if exe is None:
        raise LatexNotFoundError(
            "lualatex が見つかりません。MiKTeX か TeX Live を入れてください"
            "（Windows: winget install MiKTeX.MiKTeX）。"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    figdir_tex = figures_dir.resolve().as_posix().rstrip("/") + "/"
    # \def を先頭に置くと、main.tex の \physlabfigdir の既定値より優先される
    entry = rf"\def\physlabfigdir{{{figdir_tex}}}\input{{{main}}}"
    cmd = [exe, "-interaction=nonstopmode", "-halt-on-error", "-file-line-error"]
    if "miktex" in exe.lower():
        cmd.append("--enable-installer")  # 足りないパッケージは確認せずに自動で入れる（初回のみ）
    cmd += [f"-output-directory={out_dir.resolve()}", f"-jobname={jobname}", entry]

    log_path = out_dir / f"{jobname}.log"
    summary = summarize_log("")
    for n in range(1, max_passes + 1):
        proc = subprocess.run(  # noqa: S603 - 引数は固定の組み立てで、ユーザー入力を含まない
            cmd, cwd=tex_dir, capture_output=True, text=True, encoding="utf-8", errors="replace"
        )
        text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
        if proc.returncode != 0:
            raise LatexBuildError(f"lualatex が失敗しました（{n} 回目）:\n{error_excerpt(text)}")
        summary = summarize_log(text)
        if not summary.needs_rerun:
            return BuildResult(out_dir / f"{jobname}.pdf", log_path, n, summary)
    return BuildResult(out_dir / f"{jobname}.pdf", log_path, max_passes, summary)
