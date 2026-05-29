# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join('..', '..')))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'taster'
copyright = '2026, Luís Borges Araújo, Pablo Cardona Perez'
author = 'Luís Borges Araújo, Pablo Cardona Perez'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'numpydoc',             # NumPy documentation
    'sphinx.ext.viewcode',  # Link to local code
    'myst_nb',              # Jupyter notebooks & Markdown (superset of myst_parser)
    'sphinx_design',        # To add buttons and cards
]

templates_path = ['_templates']
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
    "private-members": True
    }

# No document TypeHints
autodoc_typehints = "none"

# Autosummary
autosummary_generate = True
autosummary_generate_overwrite = True

# MyST / myst-nb
myst_heading_anchors = 4
nb_execution_mode = "off"           # Never execute notebooks during build
nb_remove_code_outputs = True       # Strip all outputs at build time

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_show_sourcelink = False

# -- Theme configuration -----------------------------------------------------

# Sidebar configuration

html_sidebars = {
    "**": ["search-field.html", "sidebar-nav-bs.html"],
    'index': []
    }


# General theme options

html_theme_options = {
    # Logo
    'logo': {'text': project},
    # Upper bar icons
    'navbar_end': ['theme-switcher', 'navbar-icon-links'],
    # Icon links
    "icon_links": [
        # GitHub of the project
        {"name": "GitHub",
         "url": "https://github.com/Lp0lp/taster",
         "icon": "fa-brands fa-square-github",
         "type": "fontawesome",}
    ]
}

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]
