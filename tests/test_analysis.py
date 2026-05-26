import numpy as np
import pandas as pd
import pytest

from taster.analysis import TIRoutine, _format_report, process_partition


# ---------------------------------------------------------------------------
# TIRoutine — uses real XVG data from tests/data/
# ---------------------------------------------------------------------------

def test_tiroutine_mbar_cpa3_water_rep1(data_dir):
    dG, err = TIRoutine(data_dir / "CPA3" / "1" / "water", estimator="MBAR")
    assert np.isfinite(dG)
    assert np.isfinite(err)


def test_tiroutine_ti_cpa3_water_rep1(data_dir):
    dG, err = TIRoutine(data_dir / "CPA3" / "1" / "water", estimator="TI")
    assert np.isfinite(dG)
    assert np.isfinite(err)


def test_tiroutine_both_returns_nested_tuple(data_dir):
    result = TIRoutine(data_dir / "CPA3" / "1" / "water", estimator="both")
    assert len(result) == 2
    for pair in result:
        assert len(pair) == 2
        assert all(np.isfinite(v) for v in pair)


def test_tiroutine_dG_physically_plausible(data_dir):
    dG, err = TIRoutine(data_dir / "CPA3" / "1" / "water", estimator="MBAR")
    assert -200 < dG < 200
    assert err > 0


def test_tiroutine_missing_xvg_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        TIRoutine(tmp_path / "nonexistent")


# ---------------------------------------------------------------------------
# process_partition — uses symlinked tmp dirs so writes stay in tmp_path
# ---------------------------------------------------------------------------

def test_process_partition_cpa3_dataframe_shape(tmp_output_cpa3):
    df = process_partition(
        "CPA3", ["octanol-water_74-26"],
        output_dir=tmp_output_cpa3, reps=3,
    )
    assert list(df.columns) == ["rep", "solvent", "dG", "dG_err", "logP", "logP_err"]
    assert len(df) == 4  # 3 individual reps + 1 avg row


def test_process_partition_cpa3_logp_avg(tmp_output_cpa3):
    df = process_partition(
        "CPA3", ["octanol-water_74-26"],
        output_dir=tmp_output_cpa3, reps=3,
    )
    avg_logP = df[df["rep"] == "avg"]["logP"].iloc[0]
    assert avg_logP == pytest.approx(3.046, abs=0.05)


def test_process_partition_qu36_logp_avg(tmp_output_qu36):
    df = process_partition(
        "QU36", ["octanol-water_74-26"],
        output_dir=tmp_output_qu36, reps=3,
    )
    avg_logP = df[df["rep"] == "avg"]["logP"].iloc[0]
    assert avg_logP == pytest.approx(2.613, abs=0.05)


def test_process_partition_writes_results_txt(tmp_output_cpa3):
    process_partition(
        "CPA3", ["octanol-water_74-26"],
        output_dir=tmp_output_cpa3, reps=3,
    )
    assert (tmp_output_cpa3 / "CPA3" / "results.txt").exists()


def test_process_partition_npy_shape(tmp_output_cpa3):
    process_partition(
        "CPA3", ["octanol-water_74-26"],
        output_dir=tmp_output_cpa3, reps=3,
    )
    arr = np.load(tmp_output_cpa3 / "CPA3" / "partition_avg.npy")
    assert arr.shape == (1, 4)  # 1 solvent, columns: dG dG_err logP logP_err


# ---------------------------------------------------------------------------
# _format_report
# ---------------------------------------------------------------------------

def test_format_report_content():
    df = pd.DataFrame([
        {
            "rep": "1",
            "solvent": "octanol-water_74-26",
            "dG": 17.42, "dG_err": 0.055,
            "logP": 3.05, "logP_err": 0.010,
        },
        {
            "rep": "avg",
            "solvent": "octanol-water_74-26",
            "dG": 17.42, "dG_err": 0.032,
            "logP": 3.05, "logP_err": 0.006,
        },
    ])
    report = _format_report(df, "CPA3")
    assert "Molecule: CPA3" in report
    assert "Rep #1" in report
    assert "Average" in report
    assert "LogP units" in report
    assert "kJ/mol" in report
