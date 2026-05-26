import shutil
import pytest
from pathlib import Path


DATA_DIR     = Path(__file__).parent / "data"
FIXTURES_DIR = DATA_DIR / "fixtures"


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "gmx: mark test as requiring a working GROMACS installation"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as long-running (skipped by default)"
    )


# ---------------------------------------------------------------------------
# Path fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def data_dir():
    return DATA_DIR


@pytest.fixture(scope="session")
def fixtures_dir():
    return FIXTURES_DIR


# ---------------------------------------------------------------------------
# GROMACS availability (session-scoped; skips dependent tests when absent)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def gmx_executable():
    gmx = shutil.which("gmx")
    if gmx is None:
        pytest.skip("GROMACS (gmx) not found in PATH")
    return gmx


# ---------------------------------------------------------------------------
# Output-dir fixtures for process_partition tests
#
# process_partition reads  from  output_dir/resname/rep/solvent/state/fep.xvg
# and writes results.txt / partition_avg.npy into output_dir/resname/.
#
# Strategy: create a real output_dir/resname/ directory in tmp_path, then
# symlink each rep subdirectory to the corresponding tests/data directory so
# reads resolve through to real XVG data while writes stay in tmp_path.
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_output_cpa3(tmp_path, data_dir):
    mol_dir = tmp_path / "CPA3"
    mol_dir.mkdir()
    for rep in ("1", "2", "3"):
        (mol_dir / rep).symlink_to((data_dir / "CPA3" / rep).resolve())
    return tmp_path


@pytest.fixture
def tmp_output_qu36(tmp_path, data_dir):
    mol_dir = tmp_path / "QU36"
    mol_dir.mkdir()
    for rep in ("1", "2", "3"):
        (mol_dir / rep).symlink_to((data_dir / "QU36" / rep).resolve())
    return tmp_path
