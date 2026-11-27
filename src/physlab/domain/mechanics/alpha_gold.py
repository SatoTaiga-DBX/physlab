"""α 粒子を金の薄膜に当てる実験の、実際の大きさ。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章 §ラザフォードの公式

標的の反跳（金の原子核は α 粒子の約 49 倍の重さ）と、原子の電子による遮蔽は無視する。
長さは fm（1e−15 m）、エネルギーは MeV。数値は丸めた代表値。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

COULOMB_MEV_FM = 1.439964  # e² / (4π ε0) [MeV fm]
ALPHA_ELECTRON_MASS_RATIO = 7294.3  # m_α / m_e


@dataclass(frozen=True)
class AlphaOnGold:
    energy_mev: float = 5.5  # α 粒子の運動エネルギー
    foil_thickness_m: float = 4.0e-7  # 薄膜の厚さ
    charge_alpha: int = 2
    charge_gold: int = 79
    mass_number: int = 197
    atomic_radius_fm: float = 1.44e5  # 金の原子の半径（約 144 pm）
    number_density: float = 5.90e28  # 金の原子の数密度 [1/m³]

    @property
    def alpha_coefficient(self) -> float:
        """``U = α/r`` の係数 ``z Z e²/(4π ε0)`` [MeV fm]。"""
        return self.charge_alpha * self.charge_gold * COULOMB_MEV_FM

    @property
    def head_on_fm(self) -> float:
        """正面衝突の最接近距離 ``d = α/E``。"""
        return self.alpha_coefficient / self.energy_mev

    @property
    def nuclear_radius_fm(self) -> float:
        """原子核の半径の目安 ``1.2 A^{1/3}`` fm。"""
        return 1.2 * self.mass_number ** (1.0 / 3.0)

    @property
    def max_angle_from_electron(self) -> float:
        """静止した電子 1 個が α 粒子を曲げられる角の上限 ``arcsin(m_e/m_α)``。"""
        return float(np.arcsin(1.0 / ALPHA_ELECTRON_MASS_RATIO))

    def cross_section_beyond_m2(self, chi: float) -> float:
        """``χ`` より大きく曲がる断面積 ``π ρ(χ)²`` [m²]。"""
        rho_m = 0.5 * self.head_on_fm / np.tan(0.5 * chi) * 1e-15
        return float(np.pi * rho_m**2)

    def fraction_beyond(self, chi: float) -> float:
        """薄膜を通る粒子のうち ``χ`` より大きく曲がる割合 ``n t σ``（1 回だけ散乱する近似）。"""
        return self.number_density * self.foil_thickness_m * self.cross_section_beyond_m2(chi)
