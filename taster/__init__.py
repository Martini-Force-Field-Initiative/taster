"""
taster — Martini CG free energy partition coefficient calculations
"""

from importlib.metadata import version

__version__ = version("taster")

__all__ = [
    "DEFAULT_STATES",
    "TIRoutine",
    "available_solvents",
    "prepare_partition_setup",
    "process_partition",
    "run_partition_workflow",
    "run_partitions",
]

from . import utils
from .analysis import TIRoutine as TIRoutine
from .analysis import process_partition as process_partition
from .prepare import prepare_partition_setup as prepare_partition_setup
from .run import DEFAULT_STATES as DEFAULT_STATES
from .run import run_partitions as run_partitions
from .workflow import run_partition_workflow as run_partition_workflow


def available_solvents():
    """
    Return the set of solvent names bundled with taster.

    Returns
    -------
    set of str
        Solvent names available for use in partition calculations.
    """
    return utils._get_available_solvents()
