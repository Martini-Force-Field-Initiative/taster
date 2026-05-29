"""High-level orchestration of the full partition coefficient workflow.

This module exposes ``run_partition_workflow``, the single entry point that
sequentially calls :mod:`taster.prepare`, :mod:`taster.run`, and
:mod:`taster.analysis` to produce LogP values from a molecule ITP and CG
structure file.
"""
from .prepare import prepare_partition_setup
from .run import run_partitions
from .analysis import process_partition
from .utils import _get_available_solvents


def run_partition_workflow(itp, structure, solvents=None, reference='water',
                           T=298, reps=3, ncores=None, output_dir='./Partitions',
                           gmx='gmx', cutoff=5000, estimator='MBAR'):
    """
    Run the full partition coefficient workflow: prepare, run, and analyse.
 
    Parameters
    ----------
    itp : str or Path
        Path to the molecule ITP file.
    structure : str or Path
        Path to the input CG structure (must contain exactly one residue type).
    solvents : list of str, optional
        All solvent names including the reference (e.g. ['water', 'octanol', 'hexadecane']).
        Defaults to all available solvents in taster.data.solvents.
    reference : str, optional
        Name of the reference solvent to compute partitioning against.
        Must be present in solvents. Defaults to 'water'.
    T : float
        Temperature (K) at which the partitioning will be run/calculated.
    reps : int, optional
        Number of replicates. Defaults to 3.
    ncores : int or None, optional
        Number of parallel processes. Defaults to None (auto-detect).
    output_dir : str or Path, optional
        Root directory for all output. Defaults to './Partitions'.
    gmx : str, optional
        GROMACS executable name or path. Defaults to 'gmx'.
    cutoff : int, optional
        Number of initial frames to discard as equilibration. Defaults to 5000.
    estimator : str, optional
        Free energy estimator, 'TI' or 'MBAR'. Defaults to 'MBAR'.
 
    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: rep, solvent, dG, dG_err, logP, logP_err.
    """
    available = _get_available_solvents()
 
    if solvents is None:
        solvents = sorted(available)
 
    invalid = set(solvents) - available
    if invalid:
        raise ValueError(f"Unknown solvents: {sorted(invalid)}. "
                         f"Available: {sorted(available)}")
 
    if reference not in solvents:
        raise ValueError(f"Reference solvent '{reference}' must be in solvents list: {solvents}")
 
    organic_solvents = [s for s in solvents if s != reference]
 
    resname = prepare_partition_setup(itp, structure, solvents,
                                      reps=reps, output_dir=output_dir, gmx=gmx)
 
    run_partitions(resname, solvents, reps=reps, ncores=ncores,
                   output_dir=output_dir, gmx=gmx, T=T)
 
    return process_partition(resname, organic_solvents, water=reference,
                             reps=reps, output_dir=output_dir,
                             T=T, cutoff=cutoff, estimator=estimator)