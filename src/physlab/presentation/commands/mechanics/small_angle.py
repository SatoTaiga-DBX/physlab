"""`physlab run small-angle`: 微小角散乱。直線の経路に沿った力積の近似を確かめる。

設計書: docs/design/007-collisions.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.small_angle import (
    FAMILIES,
    error_sweep,
    impulse_demo,
    small_angle_tail,
)
from physlab.infrastructure.rendering.mechanics.small_angle import (
    render_error_figure,
    render_impulse_figure,
)
from physlab.infrastructure.storage import figure_path

NAME = "small-angle"
TITLE = "微小角散乱：まっすぐ進むとみなして力積を計算する"
POWERS = (1.0, 2.0, 4.0)


def figures() -> list[Path]:
    """教科書の図 2 枚。"""
    sweeps = {name: error_sweep(name) for name in FAMILIES}
    tails = [small_angle_tail(n) for n in POWERS]
    return [
        render_impulse_figure(impulse_demo(), figure_path("sma-impulse")),
        render_error_figure(sweeps, tails, figure_path("sma-error")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.parse_args(argv)

    demo = impulse_demo()
    print(f"[{NAME}] ガウス型の丘を E/U0 = 10、ρ = {demo.rho:g} a で通り抜ける")
    print(
        f"  偏向角: 厳密 {np.degrees(demo.exact_angle):.4f}°、"
        f"力積の近似 {np.degrees(demo.approx_angle):.4f}°"
    )
    print("  同じ ρ（= a）で力の強さを変える（相対誤差は、角の何乗で減るか）")
    for name in FAMILIES:
        sweep = error_sweep(name)
        print(
            f"    {name:15s} θ = {sweep.exact[0]:.4f} rad で誤差 {sweep.relative_error[0]:.1e}、"
            f"θ = {sweep.exact[-1]:.3f} rad で {sweep.relative_error[-1]:.1e}、"
            f"傾き {sweep.slope:.3f}"
        )
    print("  U = α/rⁿ の小さい角の断面積（理論の傾き −2 − 2/n）")
    for n in POWERS:
        tail = small_angle_tail(n)
        print(f"    n = {n:g}: 傾き {tail.slope:.4f}（理論 {-2 - 2 / n:g}）")

    print("図を作成中...")
    for path in figures():
        print(f"  -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
