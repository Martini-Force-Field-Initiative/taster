# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join("..", "..")))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "martini-taster"
copyright = "2026, Luís Borges Araújo, Pablo Cardona Perez"
author = "Luís Borges Araújo, Pablo Cardona Perez"

# Read the project version from pyproject.toml (PEP 621 `project.version`
# or Poetry `tool.poetry.version`). Do not silently default — fail loudly so
# docs builds surface misconfiguration.
from pathlib import Path

pyproject_file = Path(__file__).resolve().parents[2] / "pyproject.toml"
if not pyproject_file.exists():
    raise FileNotFoundError(f"pyproject.toml not found at {pyproject_file}")

pyproject_text = pyproject_file.read_text(encoding="utf-8")

try:
    import tomllib as _tomllib
except ImportError:
    try:
        import toml as _tomllib
    except ImportError:
        _tomllib = None

if _tomllib is None:
    raise RuntimeError(
        "Unable to parse pyproject.toml: Use Python 3.11+ or the 'toml' package"
    )

# Determine parser-specific decode exception classes and catch them explicitly
if hasattr(_tomllib, "TOMLDecodeError"):
    _decode_exc_types = (_tomllib.TOMLDecodeError,)
elif hasattr(_tomllib, "TomlDecodeError"):
    _decode_exc_types = (_tomllib.TomlDecodeError,)
else:
    # Conservative fallback: some parsers may raise ValueError on bad input
    _decode_exc_types = (ValueError,)

try:
    data = _tomllib.loads(pyproject_text)
except _decode_exc_types as exc:
    raise RuntimeError(f"Failed to parse pyproject.toml: {exc}") from exc

version = data.get("project", {}).get("version") or data.get("tool", {}).get(
    "poetry", {}
).get("version")
if not version:
    raise RuntimeError(
        "Version not found in pyproject.toml: expected 'project.version' or 'tool.poetry.version'"
    )

release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "numpydoc",  # NumPy documentation
    "sphinx.ext.viewcode",  # Link to local code
    "myst_nb",  # Jupyter notebooks & Markdown (superset of myst_parser)
    "sphinx_design",  # To add buttons and cards
]

templates_path = ["_templates"]
exclude_patterns = []
language = "en"

# -- Extension config --------------------------------------------------------

# Numpydoc
numpydoc_show_class_members = True
numpydoc_class_members_toctree = False
numpydoc_show_inherited_class_members = False

# Autodoc
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "private-members": True,
}

# No document TypeHints
autodoc_typehints = "none"

# Autosummary
autosummary_generate = True
autosummary_generate_overwrite = True

# MyST / myst-nb
myst_heading_anchors = 4
nb_execution_mode = "off"  # Never execute notebooks during build
nb_remove_code_outputs = True  # Strip all outputs at build time

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_show_sourcelink = False

# -- Theme configuration -----------------------------------------------------

# Sidebar configuration

html_sidebars = {"**": ["search-field.html", "sidebar-nav-bs.html"], "index": []}


# General theme options

html_theme_options = {
    # Logo
    "logo": {"text": project},
    # Upper bar icons
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    # Icon links
    "icon_links": [
        # GitHub of the project
        {
            "name": "GitHub",
            "url": "https://github.com/Martini-Force-Field-Initiative/taster",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        }
    ],
}

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]
