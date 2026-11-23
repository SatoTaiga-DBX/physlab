"""CLI の入出力と、コマンドの自動検出。"""

from pathlib import Path

from physlab.presentation.cli import discover, main


def test_discover_finds_pendulum() -> None:
    sims = discover()
    assert "pendulum" in sims
    assert sims["pendulum"].TITLE


def test_list_prints_pendulum(capsys) -> None:
    assert main(["list"]) == 0
    assert "pendulum" in capsys.readouterr().out


def test_unknown_simulation_returns_error_code(capsys) -> None:
    assert main(["run", "no-such-simulation"]) == 2
    assert "pendulum" in capsys.readouterr().err  # 候補を案内する


def test_run_pendulum_writes_outputs(tmp_path, capsys) -> None:
    code = main(["run", "pendulum", "--duration", "0.5", "--fps", "10", "--energy-duration", "2"])
    assert code == 0
    written = {p.name for p in (tmp_path / "outputs" / "mechanics").iterdir()}
    assert "pendulum_energy_drift.png" in written
    assert any(Path(n).suffix in {".mp4", ".gif"} for n in written)
    assert "最大絶対誤差" in capsys.readouterr().out
