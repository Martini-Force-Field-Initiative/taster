import sys

import pytest

from taster.cli import main


def test_cli_help_exits_zero(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["taster", "--help"])
    with pytest.raises(SystemExit) as exc_info:
        main()
    assert exc_info.value.code == 0


def test_cli_missing_required_args(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["taster"])
    with pytest.raises(SystemExit) as exc_info:
        main()
    assert exc_info.value.code == 2


def test_cli_missing_structure_arg(monkeypatch, tmp_path):
    monkeypatch.setattr(sys, "argv", ["taster", "--itp", str(tmp_path / "mol.itp")])
    with pytest.raises(SystemExit) as exc_info:
        main()
    assert exc_info.value.code == 2


def test_cli_error_handling_exit_code(monkeypatch, tmp_path):
    # Pass nonexistent files: argparse succeeds, workflow raises, CLI exits 1.
    monkeypatch.setattr(sys, "argv", [
        "taster",
        "--itp", str(tmp_path / "nonexistent.itp"),
        "--structure", str(tmp_path / "nonexistent.gro"),
        "--solvents", "water",
        "--output-dir", str(tmp_path),
    ])
    with pytest.raises(SystemExit) as exc_info:
        main()
    assert exc_info.value.code == 1
