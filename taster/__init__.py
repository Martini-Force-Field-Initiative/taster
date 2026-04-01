"""
taster — Martini CG free energy partition coefficient calculations
"""
 
__version__ = "0.1.0"

from .workflow import run_partition_workflow
from .prepare import prepare_partition_setup
from .run import run_partitions
from .analysis import process_partition, TIRoutine 
 
def available_solvents():
    """
    Return the set of solvent names bundled with taster.
 
    Returns
    -------
    set of str
        Solvent names available for use in partition calculations.
    """
    return _get_available_solvents()
 