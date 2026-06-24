"""Parallel execution of Martini FEP/TI simulations via GROMACS.

This module drives the minimisation, relaxation, and FEP production runs for
every lambda state, solvent, and replicate combination using Python
multiprocessing, with a semaphore capping concurrency and a pool of CPU pin
offsets shared across all running tasks.
"""
import os
from importlib.resources import files
from pathlib import Path
from multiprocessing import Process, Semaphore, Queue

from .utils import _run, _replace_words_in_file


DEFAULT_STATES = list(range(12))


def _run_ti_state(resname, state, workingdir, offset=0, gmx='gmx', T=298,
                  fep_min_mdp=None, fep_rel_mdp=None, fep_prod_mdp=None):
    """
    Run minimization, relaxation, and FEP production for a single lambda state.

    Parameters
    ----------
    resname : str
        Residue name of the molecule.
    state : int
        Lambda state index.
    workingdir : str or Path
        Path to the solvent directory containing system.gro and system.top.
    T : float
        Temperature (K) at which simulations will be run.
    offset : int, optional
        CPU pin offset for mdrun. Defaults to 0.
    gmx : str, optional
        GROMACS executable name or path. Defaults to 'gmx'.
    fep_min_mdp, fep_rel_mdp, fep_prod_mdp : str or Path, optional
        Paths to MDP template files. Defaults to bundled taster templates.
    """
    if fep_min_mdp is None:  fep_min_mdp  = files("taster.data.mdps") / "template-CG_minimization-fep.mdp"
    if fep_rel_mdp is None:  fep_rel_mdp  = files("taster.data.mdps") / "template-CG_relaxation-fep.mdp"
    if fep_prod_mdp is None: fep_prod_mdp = files("taster.data.mdps") / "template-CG_fep-vdw.mdp"

    workingdir = Path(workingdir).resolve()
    state_dir  = workingdir / str(state)
    state_dir.mkdir(exist_ok=True)

    # Write MDP files with molecule name and lambda state substituted
    _replace_words_in_file(fep_min_mdp,  state_dir / 'min.mdp', ['MOL', 'INIT-LAMBDA-STATE'], [resname, str(state)])
    _replace_words_in_file(fep_rel_mdp,  state_dir / 'rel.mdp', ['MOL', 'INIT-LAMBDA-STATE', 'TEMPERATURE'], [resname, str(state), str(T)])
    _replace_words_in_file(fep_prod_mdp, state_dir / 'fep.mdp', ['MOL', 'INIT-LAMBDA-STATE', 'TEMPERATURE'], [resname, str(state), str(T)])

    system_gro = workingdir / 'system.gro'
    system_top = workingdir / 'system.top'

    with open(state_dir / 'log.out', 'w') as log:
        try:
            # Minimization
            _run([gmx, 'grompp', '-f', 'min.mdp', '-c', str(system_gro),
                  '-p', str(system_top), '-o', 'min.tpr', '-maxwarn', '2'],
                 log=log, cwd=state_dir)
            _run([gmx, 'mdrun', '-deffnm', 'min', '-v', '-nt', '1',
                  '-pin', 'on', '-pinoffset', str(offset),
                  '-pme', 'cpu', '-pmefft', 'cpu', '-bonded', 'cpu',
                  '-update', 'cpu', '-nb', 'cpu'],
                 log=log, cwd=state_dir)

            # Relaxation
            _run([gmx, 'grompp', '-f', 'rel.mdp', '-c', 'min.gro',
                  '-p', str(system_top), '-o', 'rel.tpr', '-maxwarn', '3'],
                 log=log, cwd=state_dir)
            _run([gmx, 'mdrun', '-deffnm', 'rel', '-v', '-nt', '1',
                  '-pin', 'on', '-pinoffset', str(offset),
                  '-pme', 'cpu', '-pmefft', 'cpu', '-bonded', 'cpu',
                  '-update', 'cpu', '-nb', 'cpu'],
                 log=log, cwd=state_dir)

            # Cleanup relaxation files
            for pattern in ['mdout.mdp', 'min.*', 'rel.cpt', 'rel.edr',
                            'rel.log', 'rel.mdp', 'rel.tpr', 'rel.xtc', 'rel.xvg']:
                for f in state_dir.glob(pattern):
                    f.unlink(missing_ok=True)

            # FEP production
            _run([gmx, 'grompp', '-f', 'fep.mdp', '-c', 'rel.gro',
                  '-p', str(system_top), '-o', 'fep.tpr', '-maxwarn', '3'],
                 log=log, cwd=state_dir)
            _run([gmx, 'mdrun', '-deffnm', 'fep', '-v', '-nt', '1',
                  '-pin', 'on', '-pinoffset', str(offset),
                  '-pme', 'cpu', '-pmefft', 'cpu', '-bonded', 'cpu',
                  '-update', 'cpu', '-nb', 'cpu'],
                 log=log, cwd=state_dir)

        except Exception as e:
            log.write(f"\n ERROR in state {state}: {e}\n")
            raise


def _tracked_ti_state(resname, state, workingdir, offset, gmx, sem, offset_pool, T=298):
    """
    Wrapper around _run_ti_state that releases the semaphore slot and CPU
    pin offset on completion.

    Parameters
    ----------
    sem : multiprocessing.Semaphore
        Semaphore to release when the state finishes.
    offset_pool : multiprocessing.Queue
        Pool of free CPU pin offsets; `offset` is returned to it when the
        state finishes, so it's only ever reused once actually free.
    """
    try:
        _run_ti_state(resname, state, workingdir, offset=offset, gmx=gmx, T=T)
    finally:
        offset_pool.put(offset)
        sem.release()


def run_partitions(resname, solvents, reps=3, T=298,
                   output_dir='./Partitions', ncores=None, gmx='gmx',
                   states=DEFAULT_STATES):
    """
    Run TI simulations for all lambda states, solvents, and replicates locally
    using multiprocessing, with a semaphore capping concurrency and a pool of
    CPU pin offsets handed out to whichever task starts next.

    Parameters
    ----------
    resname : str
        Residue name of the molecule, used to locate the partition directory.
    solvents : list of str
        Solvent names to run.
    reps : int, optional
        Number of replicates. Defaults to 3.
    T : float
        Temperature (K) at which simulations will be run.
    output_dir : str or Path, optional
        Root directory containing prepared partition files. Defaults to './Partitions'.
    ncores : int or None, optional
        Maximum number of concurrent processes. Defaults to None (auto-detect
        via os.cpu_count()).
    gmx : str, optional
        GROMACS executable name or path. Defaults to 'gmx'.
    states : list of int, optional
        Lambda states to run. Defaults to DEFAULT_STATES (0-11).
    """
    output_dir = Path(output_dir).resolve()
    if ncores is None:
        ncores = os.cpu_count() or 1
    sem         = Semaphore(ncores)
    offset_pool = Queue()
    for i in range(ncores):
        offset_pool.put(i)
    processes   = []

    for rep in range(1, reps + 1):
        for solvent in solvents:
            workingdir = output_dir / resname / str(rep) / solvent
            for state in states:
                sem.acquire()
                offset = offset_pool.get()
                proc = Process(target=_tracked_ti_state,
                               args=(resname, state, workingdir, offset, gmx, sem, offset_pool, T))
                proc.start()
                processes.append(proc)

    for proc in processes:
        proc.join()