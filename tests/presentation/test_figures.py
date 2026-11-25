"""教科書の図の生成（`physlab run figures`）と、textbook コマンドの入口の検査。"""

import inspect
from pathlib import Path

from physlab.presentation.cli import discover, main
from physlab.presentation.registry import figure_providers

EXPECTED_FIGURES = {
    "pendulum-period-ratio",
    "pendulum-phase-portrait",
    "pendulum-energy-drift",
    "gc-same-q",
    "gc-double-pendulum-config",
    "la-path-family",
    "la-second-variation",
    "la-curvature",
    "la-landscape",
    "la-relaxation",
    "gal-train",
    "gal-action-shift",
    "gal-two-body",
    "cons-symmetry",
    "cons-integrators",
    "cons-rotating-rod",
    "com-throw",
    "com-decomposition",
    "am-areal",
    "sim-periods",
    "sim-collapse",
    "sim-virial",
}


def test_new_commands_are_discovered_and_provide_figures() -> None:
    names = set(discover())
    assert {"least-action", "generalized-coordinates", "galilean", "figures", "textbook"} <= names
    assert {"symmetry", "center-of-mass", "areal-velocity", "similarity"} <= names
    assert set(figure_providers()) == {
        "pendulum",
        "least-action",
        "generalized-coordinates",
        "galilean",
        "symmetry",
        "center-of-mass",
        "areal-velocity",
        "similarity",
    }


def test_run_figures_writes_every_textbook_figure(tmp_path: Path, capsys) -> None:
    assert main(["run", "figures"]) == 0
    written = {p.stem for p in (tmp_path / "outputs" / "figures").glob("*.pdf")}
    assert written == EXPECTED_FIGURES
    assert all(p.stat().st_size > 1000 for p in (tmp_path / "outputs" / "figures").glob("*.pdf"))
    assert f"{len(EXPECTED_FIGURES)} 個の図" in capsys.readouterr().out


def test_run_figures_with_a_name_makes_only_that_command(tmp_path: Path) -> None:
    assert main(["run", "figures", "galilean"]) == 0
    written = {p.stem for p in (tmp_path / "outputs" / "figures").glob("*.pdf")}
    assert written == {"gal-train", "gal-action-shift", "gal-two-body"}


def test_figure_functions_take_no_arguments() -> None:
    """コマンドの引数を変えても、教科書の図が別の条件で上書きされないようにする。"""
    for name, module in figure_providers().items():
        assert not inspect.signature(module.figures).parameters, name


def test_run_figures_rejects_an_unknown_command(capsys) -> None:
    assert main(["run", "figures", "no-such-command"]) == 2
    assert "候補" in capsys.readouterr().out


def test_textbook_requires_main_tex(tmp_path: Path, capsys) -> None:
    code = main(["run", "textbook", "--tex-dir", str(tmp_path / "missing"), "--no-figures"])
    assert code == 2
    assert "main.tex" in capsys.readouterr().err


def test_textbook_fails_before_compiling_when_a_figure_is_missing(tmp_path: Path, capsys) -> None:
    tex = tmp_path / "tex"
    tex.mkdir()
    (tex / "main.tex").write_text("\\physfig{no-such-figure}{c}{l}\n", encoding="utf-8")
    assert main(["run", "textbook", "--tex-dir", str(tex), "--no-figures"]) == 1
    assert "no-such-figure" in capsys.readouterr().err
