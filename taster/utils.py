"""Internal utility helpers used across the taster package.

Contains thin wrappers around subprocess for running GROMACS commands,
file-manipulation helpers for MDP template generation, and functions for
locating bundled data and XVG output files.
"""
import subprocess
from pathlib import Path
from importlib.resources import files


def _run(cmd, *, log=None, env=None, input_text=None, cwd=None):
    '''
    Short helper to assist when using subprocess to run gmx.
    '''
    subprocess.run(cmd, input=input_text,
        text=True if input_text is not None else False,
        stdout=log, stderr=subprocess.STDOUT, #if log is None goes to term.
        env=env, cwd=cwd, check=True,)


def _get_available_solvents():
    """Return the set of solvent names available in taster.data.solvents."""
    return {f.name.replace('.gro', '')
            for f in files('taster.data.solvents').iterdir()
            if f.name.endswith('.gro')}


def _replace_words_in_file(original_file_path, new_file_path,
                           words_to_replace, replacement_words):
    """
    Copy a file replacing multiple substrings in one pass.

    Parameters
    ----------
    original_file_path : str or Path
        Path to the source file.
    new_file_path : str or Path
        Path to write the modified file to.
    words_to_replace : list of str
        Substrings to search for.
    replacement_words : list of str
        Corresponding replacement strings.
    """
    content = Path(original_file_path).read_text()
    for old, new in zip(words_to_replace, replacement_words, strict=True):
        content = content.replace(old, new)
    Path(new_file_path).write_text(content)


def _get_moleculetype_name(itp_path):
    """
    Parse an ITP file and return the molecule name from its [ moleculetype ] section.

    Parameters
    ----------
    itp_path : str or Path
        Path to the ITP file.

    Returns
    -------
    str
        The molecule name (first field of the first data line under
        [ moleculetype ]).
    """
    in_section = False
    for line in Path(itp_path).read_text().splitlines():
        stripped = line.split(';', 1)[0].strip()
        if not stripped:
            continue
        if stripped.startswith('['):
            in_section = stripped.strip('[] ').lower() == 'moleculetype'
            continue
        if in_section:
            return stripped.split()[0]
    raise ValueError(f"No [ moleculetype ] section found in {itp_path}")