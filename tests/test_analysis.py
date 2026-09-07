"""Tests for taster.analysis. The TIRoutine integration tests use alchemtest's
bundled real GROMACS-format benzene dataset rather than mocking alchemlyb,
since it's tiny and fast to fit; everything else mocks TIRoutine to isolate
taster's own arithmetic from alchemlyb."""
import bz2
import shutil
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from taster.analysis import _format_report, process_partition, TIRoutine, RT, LN10


def test_format_report_includes_molecule_name_and_values():
    df = pd.DataFrame([
        dict(rep="1", solvent="octanol", dG=-12.3, dG_err=0.4, logP=2.15, logP_err=0.07),
        dict(rep="avg", solvent="octanol", dG=-12.3, dG_err=0.4, logP=2.15, logP_err=0.07),
    ])
    report = _format_report(df, "MOL")
    assert "Molecule: MOL" in report
    assert "octanol/water: 2.1500 +- 0.0700 LogP units" in report
    assert "octanol/water: -12.3000 +- 0.4000 kJ/mol" in report


def test_process_partition_partition_and_average_arithmetic(tmp_path):
    # dG/error per (rep, solvent), chosen so the expected partition
    # coefficients and averages are easy to hand-check.
    fake_results = {
        ("1", "water"): (-10.0, 0.1),
        ("1", "octanol"): (-15.0, 0.2),
        ("2", "water"): (-10.0, 0.1),
        ("2", "octanol"): (-14.0, 0.2),
    }

    def fake_ti_routine(workingdir, **kwargs):
        rep = workingdir.parent.name
        solvent = workingdir.name
        return fake_results[(rep, solvent)]

    output_dir = tmp_path / "Partitions"
    for rep in ["1", "2"]:
        for solvent in ["water", "octanol"]:
            (output_dir / "MOL" / rep / solvent).mkdir(parents=True)

    with patch("taster.analysis.TIRoutine", side_effect=fake_ti_routine):
        df = process_partition("MOL", ["octanol"], water="water", reps=2,
                               output_dir=output_dir, T=300, diagnostics=False,
                               progress=False)

    kT = RT * 300
    row1 = df[(df.rep == "1") & (df.solvent == "octanol")].iloc[0]
    expected_dG1 = -15.0 - (-10.0)
    assert row1.dG == pytest.approx(expected_dG1)
    assert row1.logP == pytest.approx(expected_dG1 / (LN10 * kT))

    avg_row = df[(df.rep == "avg") & (df.solvent == "octanol")].iloc[0]
    expected_avg_dG = np.mean([-15.0 - (-10.0), -14.0 - (-10.0)])
    assert avg_row.dG == pytest.approx(expected_avg_dG)

    assert (output_dir / "MOL" / "results.txt").is_file()
    assert (output_dir / "MOL" / "partition_avg.npy").is_file()


def _prepare_benzene_legs(base_dir):
    """Copy alchemtest's bundled benzene VDW dhdl.xvg fixtures into base_dir
    for two solvent legs, returning the list of lambda state indices."""
    alchemtest_gmx = pytest.importorskip("alchemtest.gmx")
    data = alchemtest_gmx.load_benzene()
    xvg_files = data.data["VDW"]

    for solvent in ("water", "octanol"):
        for i, f in enumerate(xvg_files):
            state_dir = base_dir / solvent / str(i)
            state_dir.mkdir(parents=True)
            with bz2.open(f, "rb") as fin, open(state_dir / "fep.xvg", "wb") as fout:
                shutil.copyfileobj(fin, fout)
    return list(range(len(xvg_files)))


def test_tiroutine_mbar_and_ti_agree_on_real_data(tmp_path):
    states = _prepare_benzene_legs(tmp_path)

    dG_mbar, err_mbar = TIRoutine(tmp_path / "water", T=300, cutoff=0,
                                   states=states, estimator="MBAR")
    dG_ti, err_ti = TIRoutine(tmp_path / "water", T=300, cutoff=0,
                              states=states, estimator="TI")

    assert dG_mbar == pytest.approx(dG_ti, abs=2.0)
    assert err_mbar > 0
    assert err_ti > 0


def test_tiroutine_diagnostics_saves_expected_figures(tmp_path):
    states = _prepare_benzene_legs(tmp_path)
    diag_dir = tmp_path / "diagnostics"

    TIRoutine(tmp_path / "water", T=300, cutoff=0, states=states,
             estimator="MBAR", diagnostics_dir=diag_dir)

    assert (diag_dir / "mbar_convergence.png").is_file()
    assert (diag_dir / "mbar_overlap_matrix.png").is_file()
    assert (diag_dir / "convergence.log").is_file()


def test_process_partition_end_to_end_with_real_data(tmp_path):
    states = _prepare_benzene_legs(tmp_path / "Partitions" / "MOL" / "1")

    df = process_partition("MOL", ["octanol"], water="water", reps=1,
                           output_dir=tmp_path / "Partitions", T=300,
                           cutoff=0, states=states, estimator="MBAR",
                           diagnostics=False, progress=False)

    assert set(df["rep"]) == {"1", "avg"}
    assert np.isfinite(df["logP"]).all()
