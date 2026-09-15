"""Tests for taster.workflow's validation and parameter wiring. The three
stage functions (prepare_partition_setup, run_partitions, process_partition)
are mocked out, since their own behavior is covered in their own test modules."""

from unittest.mock import patch

import pytest

from taster.workflow import run_partition_workflow


def test_rejects_unknown_solvent():
    with pytest.raises(ValueError, match="Unknown solvents"):
        run_partition_workflow("mol.itp", "mol.gro", solvents=["not-a-real-solvent"])


def test_rejects_reference_not_in_solvents():
    with pytest.raises(ValueError, match="must be in solvents list"):
        run_partition_workflow(
            "mol.itp",
            "mol.gro",
            solvents=["water", "hexadecane"],
            reference="octanol-water_74-26",
        )


def test_progress_flag_propagates_to_run_partitions_and_process_partition():
    with (
        patch("taster.workflow.prepare_partition_setup", return_value="MOL"),
        patch("taster.workflow.run_partitions") as mock_run,
        patch("taster.workflow.process_partition") as mock_process,
    ):
        run_partition_workflow(
            "mol.itp",
            "mol.gro",
            solvents=["water", "hexadecane"],
            reference="water",
            progress=False,
        )

    assert mock_run.call_args.kwargs["progress"] is False
    assert mock_process.call_args.kwargs["progress"] is False


def test_organic_solvents_exclude_the_reference():
    with (
        patch("taster.workflow.prepare_partition_setup", return_value="MOL"),
        patch("taster.workflow.run_partitions"),
        patch("taster.workflow.process_partition") as mock_process,
    ):
        run_partition_workflow(
            "mol.itp",
            "mol.gro",
            solvents=["water", "hexadecane", "chloroform"],
            reference="water",
        )

    organic_solvents = mock_process.call_args.args[1]
    assert set(organic_solvents) == {"hexadecane", "chloroform"}
