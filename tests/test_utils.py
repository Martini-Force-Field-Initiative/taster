"""Tests for taster.utils."""

from importlib.resources import files
from pathlib import Path

import pytest

from taster.utils import (
    _default_ff_itps,
    _get_available_solvents,
    _get_moleculetype_name,
    _get_moleculetype_names,
    _replace_words_in_file,
)


def test_get_moleculetype_name_single(tmp_path):
    itp = tmp_path / "mol.itp"
    itp.write_text("[ moleculetype ]\n; molname  nrexcl\n  MOL        1\n\n[ atoms ]\n")
    assert _get_moleculetype_name(itp) == "MOL"


def test_get_moleculetype_name_ignores_trailing_comment_on_header(tmp_path):
    itp = tmp_path / "mol.itp"
    itp.write_text(
        "[ moleculetype ]   ;; some trailing comment\n"
        "; molname  nrexcl\n"
        "  PPN        1\n"
    )
    assert _get_moleculetype_name(itp) == "PPN"


def test_get_moleculetype_name_raises_if_section_missing(tmp_path):
    itp = tmp_path / "mol.itp"
    itp.write_text("[ atoms ]\n; nothing relevant here\n")
    with pytest.raises(ValueError, match=r"No \[ moleculetype \] section"):
        _get_moleculetype_name(itp)


def test_get_moleculetype_names_returns_all_in_file_order(tmp_path):
    itp = tmp_path / "multi.itp"
    itp.write_text(
        "[moleculetype]\n; molname nrexcl\nAAA 1\n\n"
        "[moleculetype]\n; molname nrexcl\nBBB 1\n\n"
    )
    assert _get_moleculetype_names(itp) == ["AAA", "BBB"]


def test_get_moleculetype_names_finds_real_bundled_collision():
    solv_itp = files("taster.data.itps") / "martini_v3.0.0_solvents_v1.itp"
    names = _get_moleculetype_names(solv_itp)
    assert "PPN" in names
    assert "W" in names


def test_replace_words_in_file(tmp_path):
    src = tmp_path / "template.txt"
    dst = tmp_path / "out.txt"
    src.write_text("Hello MOL, lambda is INIT-LAMBDA-STATE.")
    _replace_words_in_file(src, dst, ["MOL", "INIT-LAMBDA-STATE"], ["WAT", "3"])
    assert dst.read_text() == "Hello WAT, lambda is 3."


def test_get_available_solvents_matches_bundled_data():
    solvents = _get_available_solvents()
    assert {"water", "hexadecane", "chloroform", "octanol-water_74-26"} <= solvents


def test_default_ff_itps_point_to_real_bundled_files():
    ff, solv, ions = _default_ff_itps()
    assert Path(ff).is_file()
    assert Path(solv).is_file()
    assert Path(ions).is_file()
