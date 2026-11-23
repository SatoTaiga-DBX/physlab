"""全テスト共通の設定。"""

import pytest


@pytest.fixture(autouse=True)
def isolated_output_dir(tmp_path, monkeypatch):
    """出力先を一時ディレクトリに差し替え、テストが outputs/ を汚さないようにする。"""
    monkeypatch.setenv("PHYSLAB_OUTPUT_DIR", str(tmp_path / "outputs"))
