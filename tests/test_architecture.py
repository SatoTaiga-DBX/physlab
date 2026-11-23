"""層の依存ルール（docs/design/003-layered-architecture.md）が守られているかを検査する。

契約の中身は pyproject.toml の [tool.importlinter] にある。違反すると、どの import が
どの契約を破ったかが失敗メッセージに出る。
"""

from pathlib import Path

from importlinter.cli import lint_imports

PYPROJECT = Path(__file__).resolve().parents[1] / "pyproject.toml"


def test_import_contracts_are_kept() -> None:
    assert lint_imports(config_filename=str(PYPROJECT)) == 0
