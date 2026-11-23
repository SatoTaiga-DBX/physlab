"""コマンドの自動検出。``physlab.presentation.commands`` 配下のモジュールを集める。

コマンドモジュールの規約は ``physlab.presentation.commands`` を参照。
"""

from __future__ import annotations

import importlib
import pkgutil
from types import ModuleType

from physlab.presentation import commands


def discover() -> dict[str, ModuleType]:
    """commands 配下から規約（``NAME`` / ``TITLE`` / ``main``）を満たすモジュールを集める。"""
    found: dict[str, ModuleType] = {}
    for info in pkgutil.walk_packages(commands.__path__, prefix=f"{commands.__name__}."):
        if info.ispkg:
            continue
        module = importlib.import_module(info.name)
        if all(hasattr(module, attr) for attr in ("NAME", "TITLE", "main")):
            found[module.NAME] = module
    return dict(sorted(found.items()))


def figure_providers() -> dict[str, ModuleType]:
    """教科書用の図を作る関数 ``figures()`` を持つコマンドだけを返す。"""
    return {name: m for name, m in discover().items() if callable(getattr(m, "figures", None))}
