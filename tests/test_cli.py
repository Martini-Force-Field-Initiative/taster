"""Tests for taster.cli's argument parsing and error handling."""
from unittest.mock import MagicMock, patch

import pytest

from taster.cli import main


def test_reps_defaults_to_one(monkeypatch):
    with patch("taster.cli.run_partition_workflow") as mock_workflow, \
         patch("sys.argv", ["taster", "--itp", "mol.itp", "--structure", "mol.gro"]):
        mock_workflow.return_value = MagicMock()
        main()
    assert mock_workflow.call_args.kwargs["reps"] == 1


def test_estimator_rejects_invalid_choice():
    with patch("sys.argv", ["taster", "--itp", "mol.itp", "--structure", "mol.gro",
                            "--estimator", "NOT_A_REAL_ESTIMATOR"]):
        with pytest.raises(SystemExit):
            main()


def test_errors_from_workflow_are_reported_and_exit_nonzero(capsys):
    with patch("taster.cli.run_partition_workflow", side_effect=RuntimeError("boom")), \
         patch("sys.argv", ["taster", "--itp", "mol.itp", "--structure", "mol.gro"]):
        with pytest.raises(SystemExit) as excinfo:
            main()

    assert excinfo.value.code == 1
    assert "Error: boom" in capsys.readouterr().err
