import pytest
from multiprocessing import Semaphore

from taster.run import _run_ti_state, _tracked_ti_state
from taster.prepare import prepare_partition_setup


# ---------------------------------------------------------------------------
# _tracked_ti_state — semaphore release guarantees
# (uses a real semaphore but a fake _run_ti_state-equivalent via monkeypatch)
# ---------------------------------------------------------------------------

def test_tracked_ti_state_releases_semaphore_on_success(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "taster.run._run_ti_state",
        lambda *args, **kwargs: None,
    )
    sem = Semaphore(1)
    sem.acquire()
    _tracked_ti_state("TST", 0, tmp_path, offset=0, gmx="gmx", sem=sem)
    # After the call the slot should be released — we can acquire again
    assert sem.acquire(block=False), "Semaphore was not released"


def test_tracked_ti_state_releases_semaphore_on_exception(tmp_path, monkeypatch):
    def _raise(*args, **kwargs):
        raise RuntimeError("simulated failure")

    monkeypatch.setattr("taster.run._run_ti_state", _raise)
    sem = Semaphore(1)
    sem.acquire()
    with pytest.raises(RuntimeError):
        _tracked_ti_state("TST", 0, tmp_path, offset=0, gmx="gmx", sem=sem)
    assert sem.acquire(block=False), "Semaphore was not released after exception"


# ---------------------------------------------------------------------------
# _run_ti_state — full single-state simulation (requires GROMACS, marked slow)
#
# Uses short test MDPs (nsteps=500) from tests/data/fixtures/ to keep runtime
# under a minute while still exercising minimisation → relaxation → FEP.
# ---------------------------------------------------------------------------

@pytest.mark.gmx
@pytest.mark.slow
def test_run_single_state(tmp_path, fixtures_dir, gmx_executable):
    # First prepare a solvated TST-in-water box with a correct system.top
    # that references bundled Martini ITP files (no hardcoded absolute paths).
    prepare_partition_setup(
        fixtures_dir / "TST.itp", fixtures_dir / "TST.gro",
        solvents=["water"], reps=1, output_dir=tmp_path, gmx=gmx_executable,
    )
    workingdir = tmp_path / "TST" / "1" / "water"

    _run_ti_state(
        "TST", state=0, workingdir=workingdir, gmx=gmx_executable,
        fep_min_mdp=fixtures_dir / "test_min.mdp",
        fep_rel_mdp=fixtures_dir / "test_rel.mdp",
        fep_prod_mdp=fixtures_dir / "test_fep.mdp",
    )

    xvg = workingdir / "0" / "fep.xvg"
    assert xvg.exists()
    assert xvg.stat().st_size > 0
