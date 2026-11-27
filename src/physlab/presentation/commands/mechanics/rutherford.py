"""`physlab run rutherford`: ラザフォードの公式を数えて確かめ、電荷を広げた場合と比べる。

設計書: docs/design/007-collisions.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.deflection import trace_orbits
from physlab.application.mechanics.rutherford import (
    ALPHA,
    HEAD_ON,
    ball_cross_section,
    ball_curve,
    energy_transfer,
    max_ball_angle,
    point_charge_cross_section,
)
from physlab.application.mechanics.scattering import solid_angles
from physlab.domain.mechanics.alpha_gold import AlphaOnGold
from physlab.domain.mechanics.potentials import coulomb
from physlab.domain.mechanics.scattering import (
    ball_graze_angle,
    rutherford_impact_parameter,
)
from physlab.infrastructure.rendering.mechanics.rutherford import (
    render_charge_size_figure,
    render_cross_section_figure,
    render_orbit_figure,
    render_rutherford_video,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "rutherford"
TITLE = "ラザフォードの公式：点電荷なら 1/sin⁴(χ/2)、広げると大角散乱が消える"
TOPIC = "mechanics"
BALL_RADII = (0.8, 1.4, 3.0)  # R / d
PAIR_RHO = 0.7


def _repulsive_orbits(
    rho: np.ndarray, start: float | np.ndarray = 10.0, duration: float | None = None
):
    return trace_orbits(coulomb(ALPHA), 1.0, 1.0, rho, start_radius=start, duration=duration)


def figures() -> list[Path]:
    """教科書の図 3 枚。"""
    half = np.array([0.1, 0.3, PAIR_RHO, 1.2, 2.0, 3.0])
    repulsive = _repulsive_orbits(np.concatenate([-half[::-1], half]))
    attractive = trace_orbits(coulomb(-ALPHA), 1.0, 1.0, [PAIR_RHO], start_radius=10.0)
    curves = {r: ball_curve(r) for r in BALL_RADII}
    counted = {r: ball_cross_section(r) for r in BALL_RADII}
    return [
        render_orbit_figure(repulsive, attractive, HEAD_ON, figure_path("rut-orbits")),
        render_cross_section_figure(
            point_charge_cross_section(),
            energy_transfer(),
            HEAD_ON,
            figure_path("rut-cross-section"),
        ),
        render_charge_size_figure(curves, counted, HEAD_ON, figure_path("rut-charge-size")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("--ball-radius", type=float, default=1.4, help="電荷を広げる半径 R / d")
    p.add_argument("--samples", type=int, default=200_000, help="打ち込む粒子の数")
    p.add_argument("--seed", type=int, default=1, help="乱数の種")
    p.add_argument("--fps", type=int, default=30, help="動画のフレームレート")
    p.add_argument("--no-video", action="store_true", help="動画を作らず、図だけ作る")
    args = p.parse_args(argv)
    if args.ball_radius <= 0:
        p.error("--ball-radius は正の値")

    counted = point_charge_cross_section(samples=args.samples, seed=args.seed)
    rho = rutherford_impact_parameter(np.minimum(counted.edges, np.pi - 1e-12), HEAD_ON)
    theory = -np.pi * np.diff(rho**2) / solid_angles(counted.edges)
    print(f"[{NAME}] 点電荷（d = 1）に {counted.samples} 個を打ち込み、20° 以上を角度ごとに数える")
    for k in (0, 4, 8, 12, 15):
        lo, hi = np.degrees(counted.edges[k]), np.degrees(counted.edges[k + 1])
        print(
            f"  {lo:6.1f}–{hi:6.1f}°  個数 {counted.counts[k]:6d}  "
            f"dσ/dΩ = {counted.measured[k]:9.4f}（公式 {theory[k]:9.4f}）"
        )
    z = np.max(np.abs(counted.measured - theory) / counted.error)
    print(f"  公式とのずれ（統計誤差の何倍か、全区間の最大）: {z:.1f}")

    ratio = args.ball_radius
    print(f"  電荷を半径 R = {ratio:g} d の球に広げる:")
    if ratio <= 1.0:
        print("    どの軌道も球に入らないので、点電荷と同じ")
    else:
        print(f"    点電荷と同じになるのは χ < {np.degrees(ball_graze_angle(ratio)):.2f}°")
    if ratio > 1.5:
        print(f"    中心を通り抜けるので、散乱角は最大 {np.degrees(max_ball_angle(ratio)):.2f}°")

    gold = AlphaOnGold()
    spread = gold.atomic_radius_fm / gold.head_on_fm
    print(
        f"  α 粒子 {gold.energy_mev:g} MeV と金: 正面衝突の最接近 d = {gold.head_on_fm:.1f} fm、"
        f"原子核の半径の目安 {gold.nuclear_radius_fm:.1f} fm"
    )
    print(
        f"    厚さ {gold.foil_thickness_m * 1e9:g} nm の箔で 90° 以上曲がる割合 "
        f"{gold.fraction_beyond(np.pi / 2):.2e}"
    )
    print(
        f"    正電荷が原子の大きさ（R = {spread:.0f} d）に広がるとしたときの最大の角 "
        f"{np.degrees(max_ball_angle(spread)):.4f}°"
    )
    print(f"    電子 1 個が曲げられる角の上限 {np.degrees(gold.max_angle_from_electron):.4f}°")

    print("図を作成中...")
    for path in figures():
        print(f"  -> {path}")
    if not args.no_video:
        print("動画を作成中...")
        half = np.array([0.1, 0.35, 0.6, 0.9, 1.3, 1.8, 2.4, 3.1, 4.0])
        beam = np.concatenate([-half[::-1], half])
        start = np.sqrt(7.0**2 + beam**2)  # 束の先頭をそろえる
        orbits = _repulsive_orbits(beam, start=start, duration=14.0)
        video = render_rutherford_video(
            orbits, HEAD_ON, output_path(TOPIC, "rutherford", ".mp4"), fps=args.fps
        )
        print(f"  -> {video}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
