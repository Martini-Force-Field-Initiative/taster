"""Tests for taster.run's orchestration logic, using a fake gmx executable
(see conftest.fake_gmx) instead of real GROMACS. ncores is always fixed at 2,
matching the GitHub Actions runner default, so these tests are deterministic
regardless of how many cores the host actually has."""

import itertools

import pytest

from taster.run import run_partitions


def _make_system_files(workingdir):
    workingdir.mkdir(parents=True, exist_ok=True)
    (workingdir / "system.gro").write_text("fake\n")
    (workingdir / "system.top").write_text("fake\n")


def test_run_partitions_all_states_succeed(tmp_path, fake_gmx):
    output_dir = tmp_path / "Partitions"
    _make_system_files(output_dir / "MOL" / "1" / "water")

    run_partitions(
        "MOL",
        ["water"],
        reps=1,
        ncores=2,
        gmx="gmx",
        output_dir=output_dir,
        states=[0, 1, 2, 3],
        progress=False,
    )

    for state in [0, 1, 2, 3]:
        state_dir = output_dir / "MOL" / "1" / "water" / str(state)
        assert (state_dir / "fep.gro").exists()
        assert "ERROR" not in (state_dir / "log.out").read_text()


def test_run_partitions_aggregates_a_single_failure(tmp_path, fake_gmx, monkeypatch):
    output_dir = tmp_path / "Partitions"
    _make_system_files(output_dir / "MOL" / "1" / "water")
    monkeypatch.setenv("FAIL_STATE", "2")

    with pytest.raises(RuntimeError, match=r"1 of 4 TI state\(s\) failed") as excinfo:
        run_partitions(
            "MOL",
            ["water"],
            reps=1,
            ncores=2,
            gmx="gmx",
            output_dir=output_dir,
            states=[0, 1, 2, 3],
            progress=False,
        )

    assert "state 2" in str(excinfo.value)
    # The other three states should still have completed normally.
    for state in [0, 1, 3]:
        state_dir = output_dir / "MOL" / "1" / "water" / str(state)
        assert (state_dir / "fep.gro").exists()
    assert "ERROR" in (output_dir / "MOL" / "1" / "water" / "2" / "log.out").read_text()


def test_run_partitions_never_reuses_a_still_busy_pin_offset(
    tmp_path, fake_gmx, monkeypatch
):
    log_path = tmp_path / "offsets.log"
    monkeypatch.setenv("GMX_FAKE_LOG", str(log_path))
    output_dir = tmp_path / "Partitions"
    _make_system_files(output_dir / "MOL" / "1" / "water")

    run_partitions(
        "MOL",
        ["water"],
        reps=1,
        ncores=2,
        gmx="gmx",
        output_dir=output_dir,
        states=list(range(8)),
        progress=False,
    )

    intervals = {}
    for line in log_path.read_text().splitlines():
        offset, start, end = line.split()
        intervals.setdefault(offset, []).append((float(start), float(end)))

    assert intervals, "expected the fake gmx script to have logged mdrun calls"
    for offset, spans in intervals.items():
        spans.sort()
        for (_, end_prev), (start_next, _) in itertools.pairwise(spans):
            assert start_next >= end_prev, (
                f"pin offset {offset} was reused before its previous holder finished"
            )
