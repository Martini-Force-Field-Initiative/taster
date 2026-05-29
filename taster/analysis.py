"""Free energy estimation and partition coefficient analysis.

This module post-processes the GROMACS XVG output files from TI/FEP simulations,
estimates free energy differences using TI or MBAR via alchemlyb, and computes
LogP values relative to a reference solvent.
"""
import numpy as np
import alchemlyb
from pathlib import Path
from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
from alchemlyb.estimators import TI, MBAR
 
from .run import DEFAULT_STATES
from .utils import _find_xvg_files


RT     = 0.008314  # kJ/mol/K
LN10   = np.log(10)

def TIRoutine(workingdir, T=298, cutoff=5000,
              states=None, estimator='MBAR'):
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
        if estimator_name == 'MBAR':
            data = alchemlyb.concat([extract_u_nk(str(f), T=T)[cutoff:] for f in xvg_files])
            result = MBAR().fit(data)
        else:
            data = alchemlyb.concat([extract_dHdl(str(f), T=T)[cutoff:] for f in xvg_files])
            result = TI().fit(data)
        dG    = result.delta_f_.loc[0.00, 1.00] * kT
        error = result.d_delta_f_.loc[0.00, 1.00] * kT
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
 
 
def process_partition(resname, solvents, water='water', reps=3,
                      output_dir='./Partitions', T=298,
                      cutoff=5000, states=None, estimator='MBAR'):
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
        Number of replicates. Defaults to 3.
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
 
    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: rep, solvent, dG, dG_err, logP, logP_err.
        Includes per-replicate rows and a final 'avg' row per solvent.
    """
    import pandas as pd
 
    output_dir = Path(output_dir).resolve()
    mol_dir    = output_dir / resname
    kT         = RT * T
    rows       = []
 
    for rep in range(1, reps + 1):
        # Compute dG for water reference
        water_dir          = mol_dir / str(rep) / water
        dG_water, err_water = TIRoutine(water_dir, T=T, cutoff=cutoff,
                                        states=states, estimator=estimator)
 
        # Compute dG for each organic solvent and partition against water
        for solvent in solvents:
            solvent_dir      = mol_dir / str(rep) / solvent
            dG_solv, err_solv = TIRoutine(solvent_dir, T=T, cutoff=cutoff,
                                          states=states, estimator=estimator)
 
            partition_dG     = dG_solv - dG_water
            partition_dG_err = np.sqrt(err_solv**2 + err_water**2)
            logP             = partition_dG / (LN10 * kT)
            logP_err         = partition_dG_err / (LN10 * kT)
 
            rows.append(dict(rep=str(rep), solvent=solvent,
                             dG=partition_dG, dG_err=partition_dG_err,
                             logP=logP, logP_err=logP_err))
 
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