"""動画・画像の書き出し。"""

from __future__ import annotations

from pathlib import Path

from matplotlib import animation
from matplotlib.figure import Figure


def save_animation(
    anim: animation.Animation,
    path: Path,
    fps: int = 30,
    dpi: int = 100,
) -> Path:
    """アニメーションを保存して、実際に書き出したパスを返す。

    ffmpeg が使えれば mp4（H.264）、無ければ GIF にフォールバックする。
    ``path`` の拡張子は書き出し形式に合わせて置き換わる。
    """
    if animation.FFMpegWriter.isAvailable():
        writer: animation.AbstractMovieWriter = animation.FFMpegWriter(
            fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"]
        )
        path = path.with_suffix(".mp4")
    else:
        writer = animation.PillowWriter(fps=fps)
        path = path.with_suffix(".gif")

    last_pct = -1

    def progress(i: int, n: int) -> None:
        nonlocal last_pct
        pct = (i + 1) * 100 // n
        if pct // 20 != last_pct // 20:
            print(f"  書き出し中... {pct}%")
        last_pct = pct

    anim.save(path, writer=writer, dpi=dpi, progress_callback=progress)
    return path


def save_figure(fig: Figure, path: Path, dpi: int = 150) -> Path:
    """静止画（PNG など）を保存する。"""
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path
