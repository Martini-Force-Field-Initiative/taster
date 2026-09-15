"""Shared pytest fixtures for the taster test suite."""

import os
import stat
import textwrap

import pytest

# A stand-in for the gmx CLI: instead of running any real simulation, it
# touches whatever output file(s) it was asked to produce. It can simulate
# a failure for one specific lambda state via $FAIL_STATE (matched against
# the basename of its working directory), and if $GMX_FAKE_LOG is set, every
# mdrun call appends "<pinoffset> <start> <end>" to it so tests can check
# for CPU-pin overlap.
#
# Sleep duration is deliberately *deterministic*, not random: even-numbered
# states sleep much longer than odd-numbered ones. With ncores=2 dispatch
# alternates even/odd offsets in lockstep with state parity, so this forces
# the odd (fast) state of each adjacent pair to finish first, reordering
# completion relative to dispatch every time -- regardless of how much fixed
# process-startup overhead the platform's multiprocessing start method adds
# (which would otherwise swamp small random jitter, e.g. under macOS/Windows
# "spawn", where every child re-imports the whole package from scratch).
FAKE_GMX_SCRIPT = textwrap.dedent("""\
    #!/bin/bash
    prev=""
    deffnm=""
    pinoffset=""
    for arg in "$@"; do
      if [ "$prev" = "-o" ]; then touch "$arg"; fi
      if [ "$prev" = "-deffnm" ]; then deffnm="$arg"; fi
      if [ "$prev" = "-pinoffset" ]; then pinoffset="$arg"; fi
      prev="$arg"
    done

    if [ -n "$deffnm" ]; then
      start=$(python3 -c 'import time; print(time.time())')
      state=$(basename "$PWD")
      if [ "$((state % 2))" = "0" ]; then
        sleep "${GMX_FAKE_SLEEP_EVEN:-1.5}"
      else
        sleep "${GMX_FAKE_SLEEP_ODD:-0.1}"
      fi
      touch "${deffnm}.gro" "${deffnm}.log" "${deffnm}.edr" "${deffnm}.cpt" "${deffnm}.xtc" "${deffnm}.xvg"
      end=$(python3 -c 'import time; print(time.time())')
      if [ -n "$GMX_FAKE_LOG" ]; then
        echo "$pinoffset $start $end" >> "$GMX_FAKE_LOG"
      fi
    fi

    if [ -n "$FAIL_STATE" ]; then
      cwd_base=$(basename "$PWD")
      if [ "$cwd_base" = "$FAIL_STATE" ]; then
        echo "FAKE GMX FAILURE for state $FAIL_STATE" >&2
        exit 1
      fi
    fi

    exit 0
    """)


@pytest.fixture
def fake_gmx(tmp_path, monkeypatch):
    """
    Install a fake `gmx` executable on PATH that mimics grompp/mdrun's CLI
    contract (touching expected output files) without running any real
    simulation, so taster.run's orchestration logic can be tested without
    GROMACS. Returns the directory the script lives in.
    """
    bin_dir = tmp_path / "fake_gmx_bin"
    bin_dir.mkdir()
    script = bin_dir / "gmx"
    script.write_text(FAKE_GMX_SCRIPT)
    script.chmod(script.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return bin_dir


def write_gro(path, resname, n_residues=1):
    """Write a minimal single-bead-per-residue GRO file for test fixtures."""
    lines = [f"Test molecule ({resname})", str(n_residues)]
    for i in range(1, n_residues + 1):
        lines.append(
            f"{i:>5}{resname:<5}{'C1':>5}{i:>5}{float(i):>8.3f}{0.0:>8.3f}{0.0:>8.3f}"
        )
    lines.append(f"{5.0:>10.5f}{5.0:>10.5f}{5.0:>10.5f}")
    path.write_text("\n".join(lines) + "\n")


def write_itp(path, molname):
    """Write a minimal single-bead ITP file for test fixtures."""
    path.write_text(
        "[ moleculetype ]\n"
        "; molname  nrexcl\n"
        f"  {molname}        1\n"
        "\n"
        "[ atoms ]\n"
        "; nr type resnr residue atom cgnr charge mass\n"
        f"   1   C1     1   {molname}    C1    1     0\n"
    )
