import pytest

from taster.prepare import _write_topology, prepare_partition_setup


# ---------------------------------------------------------------------------
# _write_topology — no GROMACS needed (pure MDAnalysis + file I/O)
# ---------------------------------------------------------------------------

def test_write_topology_creates_file(tmp_path, fixtures_dir):
    top_out = tmp_path / "system.top"
    _write_topology(fixtures_dir / "TST.itp", fixtures_dir / "TST.gro", output=top_out)
    assert top_out.exists()


def test_write_topology_contains_molecules_section(tmp_path, fixtures_dir):
    top_out = tmp_path / "system.top"
    _write_topology(fixtures_dir / "TST.itp", fixtures_dir / "TST.gro", output=top_out)
    content = top_out.read_text()
    assert "[ molecules ]" in content
    assert "TST" in content


def test_write_topology_molecule_count(tmp_path, fixtures_dir):
    top_out = tmp_path / "system.top"
    _write_topology(fixtures_dir / "TST.itp", fixtures_dir / "TST.gro", output=top_out)
    # The GRO has exactly 1 TST residue
    lines = top_out.read_text().splitlines()
    mol_lines = [l for l in lines if l.startswith("TST")]
    assert len(mol_lines) == 1
    assert mol_lines[0].split()[1] == "1"


def test_write_topology_includes_ff(tmp_path, fixtures_dir):
    top_out = tmp_path / "system.top"
    _write_topology(fixtures_dir / "TST.itp", fixtures_dir / "TST.gro", output=top_out)
    content = top_out.read_text()
    assert "martini_v3.0.0.itp" in content
    assert "martini_v3.0.0_solvents_v1.itp" in content
    assert "martini_v3.0.0_ions_v1.itp" in content


# ---------------------------------------------------------------------------
# prepare_partition_setup — ValueError when GRO has multiple residue types
# (no GROMACS needed; validation happens before any gmx call)
# ---------------------------------------------------------------------------

def test_prepare_multiple_resnames_raises(tmp_path, fixtures_dir):
    two_res_gro = tmp_path / "two_res.gro"
    two_res_gro.write_text(
        "Two-residue test\n"
        "    2\n"
        "    1TST     BB    1   0.000   0.000   0.000\n"
        "    2SOL      W    2   1.000   0.000   0.000\n"
        "   4.00000   4.00000   4.00000\n"
    )
    with pytest.raises(ValueError, match="Expected exactly 1 residue name"):
        prepare_partition_setup(
            fixtures_dir / "TST.itp", two_res_gro,
            solvents=["water"], reps=1, output_dir=tmp_path,
        )


# ---------------------------------------------------------------------------
# prepare_partition_setup — full run (requires GROMACS)
# ---------------------------------------------------------------------------

@pytest.mark.gmx
def test_prepare_creates_directory_structure(tmp_path, fixtures_dir, gmx_executable):
    prepare_partition_setup(
        fixtures_dir / "TST.itp", fixtures_dir / "TST.gro",
        solvents=["water"], reps=1, output_dir=tmp_path, gmx=gmx_executable,
    )
    assert (tmp_path / "TST" / "1" / "water" / "system.gro").exists()
    assert (tmp_path / "TST" / "1" / "water" / "system_.gro").exists()


@pytest.mark.gmx
def test_prepare_returns_resname(tmp_path, fixtures_dir, gmx_executable):
    resname = prepare_partition_setup(
        fixtures_dir / "TST.itp", fixtures_dir / "TST.gro",
        solvents=["water"], reps=1, output_dir=tmp_path, gmx=gmx_executable,
    )
    assert resname == "TST"


@pytest.mark.gmx
def test_prepare_generates_valid_topology(tmp_path, fixtures_dir, gmx_executable):
    prepare_partition_setup(
        fixtures_dir / "TST.itp", fixtures_dir / "TST.gro",
        solvents=["water"], reps=1, output_dir=tmp_path, gmx=gmx_executable,
    )
    top = (tmp_path / "TST" / "1" / "water" / "system.top").read_text()
    assert "[ molecules ]" in top
    assert "TST" in top


@pytest.mark.gmx
def test_prepare_multiple_reps(tmp_path, fixtures_dir, gmx_executable):
    prepare_partition_setup(
        fixtures_dir / "TST.itp", fixtures_dir / "TST.gro",
        solvents=["water"], reps=2, output_dir=tmp_path, gmx=gmx_executable,
    )
    for rep in ("1", "2"):
        assert (tmp_path / "TST" / rep / "water" / "system.gro").exists()
