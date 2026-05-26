import pytest

from taster.workflow import run_partition_workflow


def test_workflow_unknown_solvent_raises(tmp_path, fixtures_dir):
    with pytest.raises(ValueError, match="Unknown solvents"):
        run_partition_workflow(
            itp=fixtures_dir / "TST.itp",
            structure=fixtures_dir / "TST.gro",
            solvents=["atlantis"],
            output_dir=tmp_path,
        )


def test_workflow_reference_not_in_solvents_raises(tmp_path, fixtures_dir):
    with pytest.raises(ValueError, match="Reference solvent"):
        run_partition_workflow(
            itp=fixtures_dir / "TST.itp",
            structure=fixtures_dir / "TST.gro",
            solvents=["water"],
            reference="hexadecane",
            output_dir=tmp_path,
        )
