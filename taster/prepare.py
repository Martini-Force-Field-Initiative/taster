"""Solvated simulation box preparation for partition coefficient calculations.

This module builds and solvates the simulation boxes and generates the GROMACS
input files (topology, structure) needed to run thermodynamic integration across
all requested solvents and replicates.
"""
import warnings
import MDAnalysis as md
from importlib.resources import files
import numpy as np
from pathlib import Path
from .utils import _run, _get_moleculetype_name

# Input structures often have no box vectors, so MDAnalysis's 
# missing-unit-cell warnings here are expected noise, not something
# to act on.
warnings.filterwarnings('ignore', message='Empty box.*', category=UserWarning)
warnings.filterwarnings('ignore', message='.*missing dimension.*', category=UserWarning)


def _write_topology(itp, structure, output='system.top',
                    FFitp=None, SolvITP=None, IonsITP=None):
    """
    Write a GROMACS topology file with Martini FF includes and molecule counts.

    The molecule name written to `[ molecules ]` is taken from the ITP's
    `[ moleculetype ]` section.

    Parameters
    ----------
    itp : str or Path
        Path to the molecule ITP file.
    structure : str or Path
        Path to the GRO file to count residues from.
    output : str or Path
        Path to write the topology file to.
    FFitp, SolvITP, IonsITP : str or Path, optional
        Paths to Martini FF ITP files. Defaults to bundled taster data.
    """
    if FFitp is None:   FFitp   = files("taster.data.itps") / "martini_v3.0.0.itp"
    if SolvITP is None: SolvITP = files("taster.data.itps") / "martini_v3.0.0_solvents_v1.itp"
    if IonsITP is None: IonsITP = files("taster.data.itps") / "martini_v3.0.0_ions_v1.itp"

    header = [
        f'#include "{FFitp}"\n',
        f'#include "{SolvITP}"\n',
        f'#include "{IonsITP}"\n',
        f'#include "{itp}"\n',
        '[ system ]\n',
        'TI system built with Taster.\n',
        '\n',
        '[ molecules ]\n',
    ]

    molname  = _get_moleculetype_name(itp)
    n_residues = len(md.Universe(str(structure)).atoms.residues)
    header.append(f"{molname}    {n_residues}\n")

    with open(output, 'w') as topout:
        for line in header:
            topout.write(line)


def _build_box(cg_inputstructure, solvent, workingdir, gmx='gmx', d=2.0, neutralize=False):
    """
    Build a solvated simulation box using GROMACS editconf and solvate.

    Parameters
    ----------
    cg_inputstructure : str or Path
        Path to the input CG structure file.
    solvent : str
        Solvent name, must match a GRO file in taster.data.solvents.
    workingdir : Path
        Directory to write output files to.
    gmx : str, optional
        GROMACS executable name or path. Defaults to 'gmx'.
    d : float, optional
        Minimum distance (nm) between the solute and the box edge. Defaults to 2.0.
    neutralize : bool, optional
        Whether to neutralize the system with genion. Defaults to False.
    """
    SolventBox = files("taster.data.solvents") / f"{solvent}.gro"
    minMDP     = files("taster.data.mdps") / "min.mdp"

    with open(workingdir / 'build.log', 'w') as log:
        _run([gmx, 'editconf',
              '-f', str(cg_inputstructure),
              '-o', str(workingdir / 'system_.gro'),
              '-d', str(d),
              '-bt', 'dodecahedron'], log=log, cwd=workingdir)

        _run([gmx, 'solvate',
              '-cp', str(workingdir / 'system_.gro'),
              '-cs', str(SolventBox),
              '-o', str(workingdir / 'system.gro'),
              '-p', str(workingdir / 'system.top')], log=log, cwd=workingdir)

        if neutralize:
            _run([gmx, 'grompp',
                  '-f', str(minMDP),
                  '-c', str(workingdir / 'system.gro'),
                  '-p', str(workingdir / 'system.top'),
                  '-o', str(workingdir / 'neutralize.tpr'),
                  '-maxwarn', '100000'], log=log, cwd=workingdir)

            solvent_resname = np.unique(md.Universe(str(SolventBox)).atoms.resnames)[0]

            _run([gmx, 'genion',
                  '-s', str(workingdir / 'neutralize.tpr'),
                  '-o', str(workingdir / 'system.gro'),
                  '-p', str(workingdir / 'system.top'),
                  '-pname', 'NA', '-pq', '+1',
                  '-nname', 'CL', '-nq', '-1',
                  '-neutral'], log=log, cwd=workingdir, input_text=solvent_resname)


def prepare_partition_setup(itp, structure,
                            solvents, reps=3,
                            output_dir='./Partitions', gmx='gmx', d=2.0):
    """
    Prepare the directory structure and input files for partition TI simulations.

    Creates one solvated box and topology per solvent per replicate, under:
    output_dir / resname / rep / solvent /
    
    Parameters
    ----------
    itp : str or Path
        Path to the molecule ITP file.
    structure : str or Path
        Path to the input CG structure (must contain exactly one residue,
        i.e. a single copy of the molecule).
    solvents : list of str
        Solvent names to prepare (must match GRO files in taster.data.solvents).
    reps : int, optional
        Number of replicates to prepare. Defaults to 3.
    output_dir : str or Path, optional
        Root directory for output. Defaults to './Partitions'.
    gmx : str, optional
        GROMACS executable name or path. Defaults to 'gmx'.
    d : float, optional
        Minimum distance (nm) between the solute and the box edge. Defaults to 2.0.

    Returns
    -------
    str
        The molecule name, taken from the ITP's `[ moleculetype ]` section.
    """
    cg_itp            = Path(itp).resolve()
    cg_inputstructure = Path(structure).resolve()
    output_dir        = Path(output_dir).resolve()

    residues = md.Universe(str(cg_inputstructure)).residues
    if len(residues) != 1:
        raise ValueError(
            f"Expected exactly 1 residue (a single copy of the molecule) in "
            f"{cg_inputstructure}, found {len(residues)}: "
            f"{list(residues.resnames)}"
        )

    resname = _get_moleculetype_name(cg_itp)
    if residues.resnames[0] != resname:
        warnings.warn(
            f"Structure residue name '{residues.resnames[0]}' does not match "
            f"the ITP moleculetype name '{resname}'; using '{resname}' as "
            f"the molecule name throughout."
        )

    for rep in range(1, reps + 1):
        for solvent in solvents:
            workingdir = output_dir / resname / str(rep) / solvent
            workingdir.mkdir(parents=True, exist_ok=True)

            _write_topology(cg_itp, cg_inputstructure,
                            output=workingdir / 'system.top')
            _build_box(cg_inputstructure, solvent, workingdir,
                       gmx=gmx, d=d, neutralize=False)

    return resname