.. image:: _static/logo-with-name.png
   :alt: Taster
   :width: 480px
   :align: center

 
.. centered::
   `Licenced with LGPLv2.1 <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html>`_

**taster** is a Python package for computing partition coefficients (LogP) 
using thermodynamic integration (TI) free energy calculations of small molecules 
parameterized with the Martini 3 coarse-grained force field 

Given a molecule's CG structure and ITP file, taster automates the full workflow:
solvation box preparation across multiple solvents, parallel FEP/TI simulations via GROMACS,
free energy estimation using TI or MBAR, and LogP calculation relative to a reference solvent.


.. toctree::
   :maxdepth: 1
   :hidden:

   _userguide/userguide.rst
   _modules/modules.rst

.. grid:: 1 1 1 2
   :gutter: 2

   .. grid-item-card::
      :text-align: center
      :shadow: sm

      **User Guide**

      ^^^^^^^^^^^^^^

      The user guide covers installation, available solvents, and step-by-step
      examples for both the Python API and the command-line interface.

      ++++++++++

      .. button-ref:: _userguide/userguide
         :color: primary
         :expand:

         To the User Guide

   .. grid-item-card::
      :text-align: center
      :shadow: sm

      **Module Reference**

      ^^^^^^^^^^^^^^^^^

      The reference guide contains a detailed description of the public modules
      and their functions.

      ++++++++++

      .. button-ref:: _modules/modules
         :color: primary
         :expand:

         To the Module Reference


Features
--------

- **Solvation box preparation** — Builds solvated Martini simulation boxes for 
  the bundled solvents (water, octanol-water 74:26, hexadecane, chloroform)
- **Parallel FEP/TI simulations** — Runs all lambda states across solvents and replicates
  concurrently via GROMACS with semaphore-based CPU pinning
- **Free energy estimation** — Supports both TI and MBAR estimators via
  `alchemlyb <https://alchemlyb.readthedocs.io>`_
- **LogP calculation** — Computes partition coefficients relative to a chosen reference
  solvent with full error propagation across replicates
- **Bundled force field data** — Ships with Martini 3.0.0 ITP files, solvent boxes, and
  MDP templates; no external data files needed
- **Python API and CLI** — Use as a library or run the ``taster`` command directly


Installation
------------

Requirements:

- Python >= 3.10
- GROMACS >= 2024.3 (For molecule topologies using type 10 angles)
- NumPy >= 1.20
- pandas >= 1.0
- Matplotlib >= 3.5
- MDAnalysis >= 2.0
- alchemlyb >= 2.0
- tqdm >= 4.60
- JAX >= 0.6.2

Install `uv <https://docs.astral.sh/uv/getting-started/installation/>`_ first,
then choose the setup for your role.

For users, install taster with the dependencies needed to run and test it:

.. code-block:: console

   git clone https://github.com/Martini-Force-Field-Initiative/taster.git
   cd taster
    uv sync --group test

The project environment is created in ``.venv`` using the versions pinned in
``uv.lock``. The ``test`` dependency group provides pytest and alchemtest. Run
the test suite with:

.. code-block:: console

    uv run pytest -v tests/

For developers, install all linting, type-checking, testing, and documentation
dependencies:

.. code-block:: console

    uv sync --all-groups

