import subprocess

import pytest

from taster.utils import (
    _find_xvg_files,
    _get_available_solvents,
    _replace_words_in_file,
    _run,
)

EXPECTED_SOLVENTS = {"water", "octanol-water_74-26", "hexadecane", "chloroform"}


# ---------------------------------------------------------------------------
# _get_available_solvents
# ---------------------------------------------------------------------------

def test_get_available_solvents_returns_set():
    assert isinstance(_get_available_solvents(), set)


def test_get_available_solvents_known_solvents():
    assert EXPECTED_SOLVENTS <= _get_available_solvents()


# ---------------------------------------------------------------------------
# _replace_words_in_file
# ---------------------------------------------------------------------------

def test_replace_words_in_file_single(tmp_path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("Hello MOL world")
    _replace_words_in_file(src, dst, ["MOL"], ["CPA3"])
    assert dst.read_text() == "Hello CPA3 world"


def test_replace_words_in_file_multiple(tmp_path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("couple-moltype = MOL\nref-t = TEMPERATURE\ninit-lambda-state = INIT-LAMBDA-STATE")
    _replace_words_in_file(
        src, dst,
        ["MOL", "TEMPERATURE", "INIT-LAMBDA-STATE"],
        ["TST", "298", "5"],
    )
    content = dst.read_text()
    assert "TST" in content
    assert "298" in content
    assert "5" in content
    assert "MOL" not in content
    assert "TEMPERATURE" not in content
    assert "INIT-LAMBDA-STATE" not in content


def test_replace_words_in_file_mismatch_raises(tmp_path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("Hello")
    with pytest.raises(ValueError):
        _replace_words_in_file(src, dst, ["a", "b"], ["x"])


# ---------------------------------------------------------------------------
# _find_xvg_files
# ---------------------------------------------------------------------------

def test_find_xvg_files_finds_all(tmp_path):
    (tmp_path / "a.xvg").write_text("@ xvg data\n")
    sub = tmp_path / "subdir"
    sub.mkdir()
    (sub / "b.xvg").write_text("@ xvg data\n")
    (sub / "other.txt").write_text("not xvg")
    result = _find_xvg_files(tmp_path)
    names = {f.name for f in result}
    assert names == {"a.xvg", "b.xvg"}


def test_find_xvg_files_skips_hidden_dirs(tmp_path):
    hidden = tmp_path / ".ipynb_checkpoints"
    hidden.mkdir()
    (hidden / "hidden.xvg").write_text("@ xvg data\n")
    (tmp_path / "visible.xvg").write_text("@ xvg data\n")
    result = _find_xvg_files(tmp_path)
    names = {f.name for f in result}
    assert "hidden.xvg" not in names
    assert "visible.xvg" in names


# ---------------------------------------------------------------------------
# _run
# ---------------------------------------------------------------------------

def test_run_subprocess_success():
    _run(["echo", "ok"])


def test_run_subprocess_failure():
    with pytest.raises(subprocess.CalledProcessError):
        _run(["false"])
