"""matplotlib の共通設定（ヘッドレス描画・日本語フォント）。"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # 画面を出さずにファイルへ描画する（動画・画像の書き出し専用）

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

# 日本語が表示できるフォントの候補（先頭から順に、入っているものを使う）
_JP_FONTS = ["Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP", "Hiragino Sans", "IPAexGothic"]


def setup_matplotlib() -> None:
    """日本語フォントと見た目のデフォルトを設定する。図を作る前に呼ぶ。"""
    installed = {f.name for f in font_manager.fontManager.ttflist}
    available = [name for name in _JP_FONTS if name in installed]
    plt.rcParams.update(
        {
            "font.family": [*available, "sans-serif"],
            "axes.unicode_minus": False,  # 日本語フォントでマイナス記号が欠けるのを防ぐ
            "axes.grid": True,
            "grid.alpha": 0.3,
            "figure.dpi": 100,
        }
    )
