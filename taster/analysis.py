"""Free energy estimation and partition coefficient analysis.

This module post-processes the GROMACS XVG output files from TI/FEP simulations,
estimates free energy differences using TI or MBAR via alchemlyb, and computes
LogP values relative to a reference solvent.
"""
import warnings
from contextlib import contextmanager

import numpy as np
import pandas as pd
import alchemlyb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm
from loguru import logger as _logger
from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
from alchemlyb.estimators import TI, MBAR
from alchemlyb.convergence import forward_backward_convergence
from alchemlyb.visualisation import (
    plot_convergence, plot_mbar_overlap_matrix, plot_ti_dhdl,
)

from .run import DEFAULT_STATES


RT     = 0.008314  # kJ/mol/K
LN10   = np.log(10)

# alchemlyb logs convergence progress via loguru, and pandas/alchemlyb raise
# FutureWarnings during fitting. Both go to stderr by default and drown out
# the terminal during convergence analysis, so route them into a log file
# next to the diagnostic figures instead.
_logger.remove()
warnings.showwarning = lambda message, category, filename, lineno, file=None, line=None: \
    _logger.warning(f"{category.__name__}: {message} ({filename}:{lineno})")


@contextmanager
def _diagnostics_log(diagnostics_dir):
    """Redirect loguru/warnings output to `<diagnostics_dir>/convergence.log`."""
    if diagnostics_dir is None:
        yield
        return
    diagnostics_dir = Path(diagnostics_dir)
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    sink_id = _logger.add(diagnostics_dir / 'convergence.log', mode='a')
    try:
        yield
    finally:
        _logger.remove(sink_id)


def _save_diagnostics(diagnostics_dir, estimator_name, data_list, result, kT):
    """
    Compute and save convergence/overlap diagnostic figures for a single
    TI or MBAR free energy estimate.

    Parameters
    ----------
    diagnostics_dir : str or Path
        Directory to save the figures in.
    estimator_name : str
        Either 'TI' or 'MBAR'.
    data_list : list of pandas.DataFrame
        Per-state dHdl (TI) or u_nk (MBAR) DataFrames used to fit `result`.
    result : alchemlyb.estimators.TI or alchemlyb.estimators.MBAR
        The fitted estimator.
    kT : float
        Thermal energy in kJ/mol used to convert alchemlyb's native kT units.
    """
    diagnostics_dir = Path(diagnostics_dir)
    prefix = estimator_name.lower()

    convergence_df = forward_backward_convergence(data_list, estimator=estimator_name)
    convergence_df = convergence_df * kT
    ax = plot_convergence(convergence_df)
    ax.set_ylabel(r'$\Delta G$ (kJ/mol)')
    ymin, ymax = ax.get_ylim()
    margin = (ymax - ymin) * 0.15
    ax.set_ylim(ymin - margin, ymax + margin)
    ax.figure.savefig(diagnostics_dir / f'{prefix}_convergence.png', dpi=300, bbox_inches='tight')
    plt.close(ax.figure)

    if estimator_name == 'MBAR':
        ax = plot_mbar_overlap_matrix(result.overlap_matrix)
        ax.figure.savefig(diagnostics_dir / 'mbar_overlap_matrix.png', dpi=300, bbox_inches='tight')
        plt.close(ax.figure)
    else:
        ax = plot_ti_dhdl(result)
        ax.figure.savefig(diagnostics_dir / 'ti_dhdl.png', dpi=300, bbox_inches='tight')
        plt.close(ax.figure)


def TIRoutine(workingdir, T=298, cutoff=5000,
              states=None, estimator='MBAR', diagnostics_dir=None):
    """
    Run TI or MBAR free energy estimation on XVG output files from a single
    solvent directory.

    Parameters
    ----------
    workingdir : str or Path
        Path to the solvent directory containing lambda state subdirectories.
    T : float, optional
        Temperature in Kelvin. Defaults to 298.
    cutoff : int, optional
        Number of initial frames to discard as equilibration. Defaults to 5000.
    states : list of int, optional
        Lambda states to include. Defaults to DEFAULT_STATES (0-11).
    estimator : str, optional
        Free energy estimator to use. One of 'TI', 'MBAR', or 'both'.
        Defaults to 'MBAR'.
    diagnostics_dir : str or Path, optional
        If given, save forward/backward convergence plots there, plus an
        MBAR overlap matrix or TI dhdl plot depending on `estimator`.
        Defaults to None (no diagnostics saved).

    Returns
    -------
    tuple
        If estimator is 'TI' or 'MBAR': (dG, error) in kJ/mol.
        If estimator is 'both': ((dG_TI, err_TI), (dG_MBAR, err_MBAR)).
    """
    if states is None:
        states = DEFAULT_STATES

    workingdir = Path(workingdir).resolve()

    # Locate dhdl xvg files directly from state subdirectories
    xvg_files = [workingdir / str(state) / 'fep.xvg' for state in states]
    missing = [f for f in xvg_files if not f.exists()]
    if missing:
        raise FileNotFoundError(f"Missing XVG files: {missing}")

    kT = RT * T

    def _extract(estimator_name):
        with _diagnostics_log(diagnostics_dir):
            if estimator_name == 'MBAR':
                data_list = [extract_u_nk(str(f), T=T)[cutoff:] for f in xvg_files]
                result = MBAR().fit(alchemlyb.concat(data_list))
            else:
                data_list = [extract_dHdl(str(f), T=T)[cutoff:] for f in xvg_files]
                result = TI().fit(alchemlyb.concat(data_list))
            dG    = result.delta_f_.loc[0.00, 1.00] * kT
            error = result.d_delta_f_.loc[0.00, 1.00] * kT
            if diagnostics_dir is not None:
                _save_diagnostics(diagnostics_dir, estimator_name, data_list, result, kT)
        return dG, error

    if estimator == 'both':
        return _extract('TI'), _extract('MBAR')
    elif estimator in ('TI', 'MBAR'):
        return _extract(estimator)
    else:
        raise ValueError(f"Unknown estimator '{estimator}'. Choose 'TI', 'MBAR', or 'both'.")

 
def _format_report(df, resname):
    """
    Format a human-readable report from a partition results DataFrame.
 
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with columns: rep, solvent, dG, dG_err, logP, logP_err.
    resname : str
        Residue name of the molecule.
 
    Returns
    -------
    str
        Formatted report string.
    """
    lines = [f'Molecule: {resname}\n']
 
    for rep, group in df[df['rep'] != 'avg'].groupby('rep'):
        lines.append(f'Rep #{rep}\n')
        for _, row in group.iterrows():
            lines.append(f"  {row['solvent']}/water: {row['logP']:.4f} +- {row['logP_err']:.4f} LogP units\n")
            lines.append(f"  {row['solvent']}/water: {row['dG']:.4f} +- {row['dG_err']:.4f} kJ/mol\n")
 
    lines.append('Average\n')
    for _, row in df[df['rep'] == 'avg'].iterrows():
        lines.append(f"  {row['solvent']}/water: {row['logP']:.4f} +- {row['logP_err']:.4f} LogP units\n")
        lines.append(f"  {row['solvent']}/water: {row['dG']:.4f} +- {row['dG_err']:.4f} kJ/mol\n")
 
    return ''.join(lines)
 
 
def process_partition(resname, solvents, water='water', reps=1,
                      output_dir='./Partitions', T=298,
                      cutoff=5000, states=None, estimator='MBAR',
                      diagnostics=True, progress=True):
    """
    Process partition coefficient TI simulations and compute LogP values.
 
    For each replicate and solvent, runs free energy estimation, computes
    the partition coefficient relative to water, then averages over replicates.
    Saves results as a text report and numpy array.
 
    Parameters
    ----------
    resname : str
        Residue name of the molecule.
    solvents : list of str
        Organic solvent names (e.g. ['octanol-water_74-26', 'hexadecane']).
        Water is handled separately via the `water` parameter.
    water : str, optional
        Name of the water solvent directory. Defaults to 'water'.
    reps : int, optional
        Number of replicates. Defaults to 1.
    output_dir : str or Path, optional
        Root directory containing simulation output. Defaults to './Partitions'.
    T : float, optional
        Temperature in Kelvin. Defaults to 298.
    cutoff : int, optional
        Number of initial frames to discard as equilibration. Defaults to 5000.
    states : list of int, optional
        Lambda states to include. Defaults to DEFAULT_STATES.
    estimator : str, optional
        Free energy estimator, 'TI' or 'MBAR'. Defaults to 'MBAR'.
    diagnostics : bool, optional
        If True, save convergence diagnostics (plus an MBAR overlap matrix
        or TI dhdl plot, depending on `estimator`) for every replicate/solvent
        leg under `<mol_dir>/diagnostics/<rep>/<solvent>/`. Defaults to True.
    progress : bool, optional
        Whether to display a tqdm progress bar. Defaults to True.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: rep, solvent, dG, dG_err, logP, logP_err.
        Includes per-replicate rows and a final 'avg' row per solvent.
    """
    output_dir = Path(output_dir).resolve()
    mol_dir    = output_dir / resname
    kT         = RT * T
    rows       = []

    def _diagnostics_dir(rep, solvent_name):
        if not diagnostics:
            return None
        return mol_dir / 'diagnostics' / str(rep) / solvent_name

    total = reps * (1 + len(solvents))
    pbar  = tqdm(total=total, desc=f"{resname} analysis", disable=not progress)

    for rep in range(1, reps + 1):
        # Compute dG for water reference
        water_dir          = mol_dir / str(rep) / water
        dG_water, err_water = TIRoutine(water_dir, T=T, cutoff=cutoff,
                                        states=states, estimator=estimator,
                                        diagnostics_dir=_diagnostics_dir(rep, water))
        pbar.update(1)

        # Compute dG for each organic solvent and partition against water
        for solvent in solvents:
            solvent_dir      = mol_dir / str(rep) / solvent
            dG_solv, err_solv = TIRoutine(solvent_dir, T=T, cutoff=cutoff,
                                          states=states, estimator=estimator,
                                          diagnostics_dir=_diagnostics_dir(rep, solvent))
            pbar.update(1)

            partition_dG     = dG_solv - dG_water
            partition_dG_err = np.sqrt(err_solv**2 + err_water**2)
            logP             = partition_dG / (LN10 * kT)
            logP_err         = partition_dG_err / (LN10 * kT)

            rows.append(dict(rep=str(rep), solvent=solvent,
                             dG=partition_dG, dG_err=partition_dG_err,
                             logP=logP, logP_err=logP_err))

    pbar.close()
 
    df = pd.DataFrame(rows)
 
    # Average over replicates
    avg_rows = []
    for solvent in solvents:
        sol_df = df[df['solvent'] == solvent]
        n      = len(sol_df)
        avg_rows.append(dict(
            rep='avg', solvent=solvent,
            dG=sol_df['dG'].mean(),
            dG_err=np.sqrt((sol_df['dG_err']**2).sum()) / n,
            logP=sol_df['logP'].mean(),
            logP_err=np.sqrt((sol_df['logP_err']**2).sum()) / n,
        ))
    df = pd.concat([df, pd.DataFrame(avg_rows)], ignore_index=True)
 
    # Save results
    report = _format_report(df, resname)
    with open(mol_dir / 'results.txt', 'w') as f:
        f.write(report)
 
    avg_array = df[df['rep'] == 'avg'][['dG', 'dG_err', 'logP', 'logP_err']].to_numpy()
    np.save(mol_dir / 'partition_avg.npy', avg_array)
 
    return df