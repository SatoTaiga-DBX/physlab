"""教科書の TeX ソースの整合性（TeX を実行せずに調べられる範囲）。

TeX 本体のコンパイルと未解決参照の検査は ``physlab run textbook`` が行う。
ここでは、人が書き間違えやすい「問題と解答の対応」「ファイルの存在」「ラベルの重複」を調べる。
"""

import re
from pathlib import Path

import pytest

TEX_DIR = Path(__file__).resolve().parents[2] / "docs/textbook/tex"


def _read(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    return re.sub(r"(?m)^\s*%.*$", "", text)  # コメント行は除く


@pytest.fixture(scope="module")
def sources() -> dict[Path, str]:
    # 雛形（_template.tex）は本の一部ではないので対象外
    return {p: _read(p) for p in sorted(TEX_DIR.rglob("*.tex")) if not p.name.startswith("_")}


def test_main_inputs_exist(sources: dict[Path, str]) -> None:
    for path, text in sources.items():
        for target in re.findall(r"\\input\{([^}]+)\}", text):
            assert (TEX_DIR / f"{target}.tex").is_file(), f"{path}: \\input{{{target}}} が無い"


def test_every_exercise_has_one_solution(sources: dict[Path, str]) -> None:
    text = "\n".join(sources.values())
    exercises = re.findall(r"\\begin\{exercise\}(?:\[[^\]]*\])?\{([^}]+)\}", text)
    solutions = re.findall(r"\\begin\{solution\}\{([^}]+)\}", text)
    assert exercises, "問題が 1 つも無い（正規表現の見直しが必要）"
    assert len(set(exercises)) == len(exercises), "問題の鍵が重複している"
    assert len(set(solutions)) == len(solutions), "解答の鍵が重複している"
    assert set(exercises) - set(solutions) == set(), "解答の無い問題がある"
    assert set(solutions) - set(exercises) == set(), "問題の無い解答がある"


def test_labels_are_unique(sources: dict[Path, str]) -> None:
    labels = re.findall(r"\\label\{([^}]+)\}", "\n".join(sources.values()))
    duplicated = {x for x in labels if labels.count(x) > 1}
    assert not duplicated, f"ラベルが重複: {sorted(duplicated)}"


def test_cited_keys_are_in_bibliography(sources: dict[Path, str]) -> None:
    text = "\n".join(sources.values())
    cited = {
        k.strip() for group in re.findall(r"\\cite\{([^}]+)\}", text) for k in group.split(",")
    }
    defined = set(re.findall(r"\\bibitem\{([^}]+)\}", text))
    assert cited <= defined, f"参考文献に無い鍵: {sorted(cited - defined)}"
