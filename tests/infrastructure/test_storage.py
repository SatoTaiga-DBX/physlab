from pathlib import Path

from physlab.infrastructure.storage import output_path


def test_output_path_follows_env_var_and_creates_parent() -> None:
    path = output_path("mechanics", "demo", ".png")
    assert path.name == "demo.png"
    assert path.parent.name == "mechanics"
    assert path.parent.is_dir()
    assert "outputs" in Path(path).parts
