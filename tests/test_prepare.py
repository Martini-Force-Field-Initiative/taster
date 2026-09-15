"""Tests for taster.prepare. No real GROMACS is invoked: _build_box is
mocked out wherever it would otherwise run editconf/solvate."""

from unittest.mock import patch

import pytest
from conftest import write_gro, write_itp

from taster.prepare import _write_topology, prepare_partition_setup


def test_write_topology_uses_itp_moleculetype_name_not_structure_resname(tmp_path):
    gro = tmp_path / "mol.gro"
    itp = tmp_path / "mol.itp"
    write_gro(gro, resname="X")  # deliberately mismatched on-disk resname
    write_itp(itp, molname="MOL")
    out = tmp_path / "system.top"

    _write_topology(itp, gro, output=out)

    last_line = out.read_text().splitlines()[-1].split()
    assert last_line[0] == "MOL"
    assert last_line[1] == "1"


def test_write_topology_counts_residues(tmp_path):
    gro = tmp_path / "mol.gro"
    itp = tmp_path / "mol.itp"
    write_gro(gro, resname="MOL", n_residues=3)
    write_itp(itp, molname="MOL")
    out = tmp_path / "system.top"

    _write_topology(itp, gro, output=out)

    last_line = out.read_text().splitlines()[-1].split()
    assert last_line == ["MOL", "3"]


def test_prepare_partition_setup_rejects_multiple_residues(tmp_path):
    gro = tmp_path / "mol.gro"
    itp = tmp_path / "mol.itp"
    write_gro(gro, resname="MOL", n_residues=2)
    write_itp(itp, molname="MOL")

    with pytest.raises(ValueError, match="Expected exactly 1 residue"):
        prepare_partition_setup(
            itp, gro, solvents=["water"], reps=1, output_dir=tmp_path / "Partitions"
        )


def test_prepare_partition_setup_warns_on_resname_mismatch(tmp_path):
    gro = tmp_path / "mol.gro"
    itp = tmp_path / "mol.itp"
    write_gro(gro, resname="X")
    write_itp(itp, molname="MOL")

    with (
        patch("taster.prepare._build_box"),
        pytest.warns(UserWarning, match="does not match"),
    ):
        resname = prepare_partition_setup(
            itp, gro, solvents=["water"], reps=1, output_dir=tmp_path / "Partitions"
        )
    assert resname == "MOL"


def test_prepare_partition_setup_warns_on_bundled_solvent_name_collision(tmp_path):
    # "PPN" is acetone's moleculetype name in the bundled solvent ITP.
    gro = tmp_path / "ppn.gro"
    itp = tmp_path / "ppn.itp"
    write_gro(gro, resname="PPN")
    write_itp(itp, molname="PPN")

    with (
        patch("taster.prepare._build_box"),
        pytest.warns(UserWarning, match="collides with a moleculetype"),
    ):
        prepare_partition_setup(
            itp, gro, solvents=["water"], reps=1, output_dir=tmp_path / "Partitions"
        )


def test_prepare_partition_setup_no_warnings_when_consistent(tmp_path, recwarn):
    gro = tmp_path / "mol.gro"
    itp = tmp_path / "mol.itp"
    write_gro(gro, resname="XYZ12")
    write_itp(itp, molname="XYZ12")

    with patch("taster.prepare._build_box"):
        resname = prepare_partition_setup(
            itp, gro, solvents=["water"], reps=1, output_dir=tmp_path / "Partitions"
        )

    assert resname == "XYZ12"
    assert len(recwarn) == 0


def test_prepare_partition_setup_creates_expected_directory_layout(tmp_path):
    gro = tmp_path / "mol.gro"
    itp = tmp_path / "mol.itp"
    write_gro(gro, resname="MOL")
    write_itp(itp, molname="MOL")
    output_dir = tmp_path / "Partitions"

    with patch("taster.prepare._build_box"):
        prepare_partition_setup(
            itp, gro, solvents=["water", "hexadecane"], reps=2, output_dir=output_dir
        )

    for rep in ("1", "2"):
        for solvent in ("water", "hexadecane"):
            assert (output_dir / "MOL" / rep / solvent / "system.top").is_file()
