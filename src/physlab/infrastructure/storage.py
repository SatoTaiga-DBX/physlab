"""出力ファイルの置き場所の規約。"""

from __future__ import annotations

import os
from pathlib import Path


def output_root() -> Path:
    """出力のルート。環境変数 ``PHYSLAB_OUTPUT_DIR`` で変更でき、既定はカレントの ``outputs/``。"""
    return Path(os.environ.get("PHYSLAB_OUTPUT_DIR", "outputs"))


def output_dir(topic: str) -> Path:
    """``outputs/<topic>/`` を返す（作成する）。"""
    path = output_root() / topic
    path.mkdir(parents=True, exist_ok=True)
    return path


def output_path(topic: str, name: str, suffix: str = "") -> Path:
    """``outputs/<topic>/<name><suffix>`` を返す（親ディレクトリは作成する）。"""
    return output_dir(topic) / f"{name}{suffix}"


# 教科書（TeX）に貼る図の置き場所。TeX 側は \physlabfigdir 経由でここを読む
FIGURES_TOPIC = "figures"


def figure_path(name: str, suffix: str = ".pdf") -> Path:
    """教科書用の図 ``outputs/figures/<name><suffix>``（既定はベクタ画像の PDF）。"""
    return output_path(FIGURES_TOPIC, name, suffix)
