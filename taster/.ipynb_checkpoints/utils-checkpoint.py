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

    
def _find_xvg_files(workingdir):
    """
    Recursively find all XVG files in a directory, ignoring hidden directories.
 
    Parameters
    ----------
    workingdir : str or Path
        Directory to search.
 
    Returns
    -------
    list of Path
        All XVG files found.
    """
    workingdir = Path(workingdir)
    return [
        f for f in workingdir.rglob('*.xvg')
        if not any(part.startswith('.') for part in f.parts)
    ]