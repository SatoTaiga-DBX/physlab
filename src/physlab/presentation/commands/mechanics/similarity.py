"""`physlab run similarity`: 力学的相似（周期の冪、相似な軌道、ビリアル定理）。

設計書: docs/design/006-conservation-laws.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

from physlab.application.mechanics.similarity import (
    collapse,
    kepler_periods,
    kepler_virial,
    sweep_periods,
    virial,
)
from physlab.infrastructure.rendering.mechanics.similarity import (
    render_collapse_figure,
    render_period_figure,
    render_virial_figure,
)
from physlab.infrastructure.storage import figure_path

NAME = "similarity"
TITLE = "力学的相似：ポテンシャルの形だけで、周期の振幅依存が決まる"


def figures() -> list[Path]:
    """教科書の図 3 枚（``outputs/figures/sim-*.pdf``）。"""
    return [
        render_period_figure(sweep_periods(), kepler_periods(), figure_path("sim-periods")),
        render_collapse_figure(collapse(), figure_path("sim-collapse")),
        render_virial_figure(virial(), figure_path("sim-virial")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.parse_args(argv)

    sweep = sweep_periods()
    print(f"[{NAME}] U = |x|^k / k の周期と振幅（A = 0.25 – 4）")
    for k, slope in sweep.slopes.items():
        error = max(abs(sweep.measured[k] / sweep.exact[k] - 1))
        theory = 1 - k / 2
        print(f"  k={k:g}: 傾き {slope:+.6f}（理論 {theory:+g}）、公式との差の最大 {error:.1e}")
    kepler = kepler_periods()
    print(f"  万有引力: 長半径に対する傾き {kepler.slope:.6f}（理論 1.5）")
    print(f"相似な軌道の重なり（k=4）: 最大の差 {collapse().max_difference:.1e}")
    ratios = ", ".join(f"k={k:g}: {v:.6f}" for k, v in virial().over_whole_periods.items())
    print(f"ビリアル比 2⟨T⟩/⟨U⟩（整数周期）: {ratios}、万有引力 {kepler_virial():.6f}")

    print("図を作成中...")
    for path in figures():
        print(f"  -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
