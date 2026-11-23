"""出力ファイルの置き場所の規約。"""

from __future__ import annotations

import os
from pathlib import Path


def output_path(topic: str, name: str, suffix: str = "") -> Path:
    """``outputs/<topic>/<name><suffix>`` を返す（親ディレクトリは作成する）。

    ルートは環境変数 ``PHYSLAB_OUTPUT_DIR`` で変更でき、既定はカレントの ``outputs/``。
    """
    root = Path(os.environ.get("PHYSLAB_OUTPUT_DIR", "outputs"))
    path = root / topic / f"{name}{suffix}"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
