"""`physlab run scattering`: 衝突パラメータから散乱角、そして断面積へ。虹散乱。

設計書: docs/design/007-collisions.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.deflection import (
    deflection,
    deflection_in_chunks,
    orbit_deflection,
    trace_orbits,
)
from physlab.application.mechanics.scattering import (
    DeflectionCurve,
    count_cross_section,
    curve_cross_section,
    deflection_curve,
    hard_sphere_angle_from_reflection,
    hard_sphere_lab_cross_section,
    hill,
    rainbow,
    rainbow_edges,
)
from physlab.domain.mechanics.scattering import hard_sphere_angle
from physlab.infrastructure.rendering.mechanics.scattering import (
    CurveSeries,
    render_beam_video,
    render_cross_section_figure,
    render_deflection_figure,
    render_lab_figure,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "scattering"
TITLE = "散乱と断面積：打ち込んで、曲がった向きを数える（虹散乱）"
TOPIC = "mechanics"
ENERGY_RATIO = 1.5  # E / U0（丘の高さ U0 = 1）
RHO_MAX = 3.5
ORBIT_RHO = np.array([0.15, 0.4, 0.8, 1.3, 1.9])


def _speed(energy: float) -> float:
    return float(np.sqrt(2.0 * energy))


def _series(label: str, energy: float, color: str) -> CurveSeries:
    curve = deflection_curve(hill(), energy, 2.6, n=800)
    orbits = trace_orbits(hill(), 1.0, _speed(energy), ORBIT_RHO, start_radius=8.0)
    return CurveSeries(label, curve, color, ORBIT_RHO, np.abs(orbit_deflection(orbits)))


def _count_hill(energy: float, samples: int, seed: int):
    curve = deflection_curve(hill(), energy, RHO_MAX)
    _, chi_r = rainbow(curve)
    edges = rainbow_edges(chi_r)
    counted = count_cross_section(
        lambda rho: deflection_in_chunks(rho, hill(), energy),
        RHO_MAX,
        edges,
        samples=samples,
        seed=seed,
    )
    return curve, chi_r, counted, curve_cross_section(curve, edges)


def _count_sphere(samples: int, seed: int):
    return count_cross_section(
        lambda rho: hard_sphere_angle_from_reflection(rho, 1.0),
        1.0,
        np.linspace(0.05, np.pi, 19),
        samples=samples,
        seed=seed,
    )


def figures() -> list[Path]:
    """教科書の図 3 枚。"""
    beam = trace_orbits(
        hill(), 1.0, _speed(ENERGY_RATIO), np.linspace(-2.5, 2.5, 41), start_radius=6.0
    )
    rho = np.linspace(1e-3, 2.6, 800)
    series = [
        CurveSeries("剛体球", DeflectionCurve(rho, hard_sphere_angle(rho, 1.0)), "0.5"),
        _series("丘 $E/U_0=0.8$", 0.8, "tab:red"),
        _series("丘 $E/U_0=1.5$", ENERGY_RATIO, "tab:blue"),
        _series("丘 $E/U_0=3$", 3.0, "tab:green"),
    ]
    curve, chi_r, counted, theory = _count_hill(ENERGY_RATIO, 200_000, 1)
    lab = [hard_sphere_lab_cross_section(m1, m2) for m1, m2 in ((1, 5), (1, 1), (5, 1))]
    return [
        render_deflection_figure(beam, series, figure_path("sca-deflection"), ENERGY_RATIO),
        render_cross_section_figure(
            _count_sphere(200_000, 1),
            counted,
            curve,
            theory,
            chi_r,
            figure_path("sca-cross-section"),
        ),
        render_lab_figure(lab, figure_path("sca-lab-frame")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("--energy", type=float, default=ENERGY_RATIO, help="E / U0（丘の高さが 1）")
    p.add_argument("--samples", type=int, default=200_000, help="打ち込む粒子の数")
    p.add_argument("--seed", type=int, default=1, help="乱数の種")
    p.add_argument("--fps", type=int, default=30, help="動画のフレームレート")
    p.add_argument("--no-video", action="store_true", help="動画を作らず、図だけ作る")
    args = p.parse_args(argv)
    if args.energy <= 0 or abs(args.energy - 1.0) < 0.05:
        p.error("--energy は正で、丘の頂上（1）から 0.05 以上離した値（頂上では粒子が止まる）")

    print(f"[{NAME}] ガウス型の丘 U = exp(−r²)、E/U0 = {args.energy:g}")
    sample_rho = np.array([0.1, 0.4, 0.8, 1.3, 2.0])
    angles = np.degrees(deflection(sample_rho, hill(), args.energy))
    for rho, angle in zip(sample_rho, angles, strict=True):
        print(f"  ρ = {rho:g} → χ = {angle:7.2f}°")
    if args.energy > 1.0:
        curve, chi_r, counted, theory = _count_hill(args.energy, args.samples, args.seed)
        rho_r, _ = rainbow(curve)
        estimate = np.sqrt(np.pi / (2.0 * np.e)) / args.energy
        used = theory > 0
        z = np.max(np.abs(counted.measured - theory)[used] / counted.error[used])
        print(
            f"  虹の角 χ_r = {np.degrees(chi_r):.2f}°（ρ = {rho_r:.3f}）、"
            f"力積の近似では {np.degrees(estimate):.2f}°"
        )
        print(f"  数えた断面積と理論のずれ（統計誤差の何倍か、最大）: {z:.1f}")
    else:
        print("  E < U0: 正面から来た粒子は跳ね返る（χ は π から単調に減り、虹の角は無い）")
    sphere = _count_sphere(args.samples, args.seed)
    print(
        f"  剛体球（a = 1）の dσ/dΩ: {sphere.measured.min():.4f} – {sphere.measured.max():.4f}"
        "（理論 0.25）"
    )

    print("図を作成中...")
    for path in figures():
        print(f"  -> {path}")
    if not args.no_video:
        print("動画を作成中...")
        speed = _speed(args.energy)
        rho = np.linspace(-3.0, 3.0, 40)
        orbits = trace_orbits(
            hill(), 1.0, speed, rho, start_radius=np.sqrt(7.0**2 + rho**2), duration=14.0 / speed
        )
        chi = np.abs(deflection(np.abs(rho), hill(), args.energy))
        curve = deflection_curve(hill(), args.energy, 3.0, n=800)
        video = render_beam_video(
            orbits, chi, curve, output_path(TOPIC, "scattering_beam", ".mp4"), fps=args.fps
        )
        print(f"  -> {video}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
