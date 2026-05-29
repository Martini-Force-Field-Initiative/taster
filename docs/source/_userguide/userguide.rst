User Guide
==========

This guide covers installation, available solvents, and step-by-step usage examples
for both the Python API and the command-line interface.


Available Solvents
------------------

taster ships with four bundled solvent boxes compatible with Martini 3:

+------------------------------+---------------------------------------------+
| Solvent name                 | Description                                 |
+==============================+=============================================+
| ``water``                    | Pure water (typical reference solvent)      |
+------------------------------+---------------------------------------------+
| ``octanol-water_74-26``      | 1-octanol / water mixture                   |
+------------------------------+---------------------------------------------+
| ``hexadecane``               | Pure hexadecane                             |
+------------------------------+---------------------------------------------+
| ``chloroform``               | Pure chloroform                             |
+------------------------------+---------------------------------------------+

You can query available solvents at runtime:

.. code-block:: python

   import taster
   print(taster.available_solvents())


Python API
----------

Full workflow example
~~~~~~~~~~~~~~~~~~~~~

The simplest entry point runs prepare, simulate, and analyse in one call:

.. code-block:: python

   import taster

   df = taster.run_partition_workflow(
       itp='MOL.itp',
       structure='MOL.gro',
       solvents=['water', 'octanol-water_74-26', 'hexadecane'],
       reference='water',
       reps=3,
       ncores=36,
       estimator='MBAR',
   )

   print(df)

The returned ``pandas.DataFrame`` has one row per (replicate, solvent) combination plus
averaged rows, with columns ``rep``, ``solvent``, ``dG``, ``dG_err``, ``logP``,
``logP_err``.


Step-by-step workflow
~~~~~~~~~~~~~~~~~~~~~

Each stage can also be called independently, which is useful when you want to inspect
intermediate results or restart from a specific step:

.. code-block:: python

   import taster

   # Step 1 — prepare solvated simulation boxes
   resname = taster.prepare_partition_setup(
       itp='MOL.itp',
       structure='MOL.gro',
       solvents=['water', 'octanol-water_74-26', 'hexadecane'],
       reps=3,
       output_dir='./Partitions',
   )

   # Step 2 — run TI simulations in parallel
   taster.run_partitions(
       resname,
       solvents=['water', 'octanol-water_74-26', 'hexadecane'],
       reps=3,
       ncores=36,
   )

   # Step 3 — analyse and compute LogP
   df = taster.process_partition(
       resname,
       solvents=['octanol-water_74-26', 'hexadecane'],
       water='water',
       reps=3,
   )

   print(df)

Command-Line Interface
----------------------

The ``taster`` command exposes the full workflow without writing any Python:

.. code-block:: bash

   taster --itp MOL.itp --structure MOL.gro \
          --solvents water octanol-water_74-26 hexadecane \
          --reference water \
          --reps 3 \
          --ncores 36 \
          --estimator MBAR

All CLI options:

.. code-block:: text

   usage: taster [-h] --itp ITP --structure STRUCTURE
                 [--solvents SOLVENTS [SOLVENTS ...]]
                 [--reference REFERENCE]
                 [--reps REPS] [--ncores NCORES]
                 [--output-dir OUTPUT_DIR] [--gmx GMX]
                 [--temperature TEMPERATURE]
                 [--cutoff CUTOFF]
                 [--estimator {TI,MBAR}]

   Compute Martini 3 partition coefficients via TI/MBAR free energy calculations.

   required arguments:
     --itp ITP             Path to the molecule ITP file
     --structure STRUCTURE Path to the input CG structure (.gro)

   optional arguments:
     --solvents            Solvent names (default: all available)
     --reference           Reference solvent for LogP (default: water)
     --reps                Number of replicates (default: 3)
     --ncores              Parallel processes (default: auto-detect)
     --output-dir          Output root directory (default: ./Partitions)
     --gmx                 GROMACS executable (default: gmx)
     --temperature         Temperature in K (default: 298)
     --cutoff              Equilibration frames to discard (default: 5000)
     --estimator           Free energy estimator: TI or MBAR (default: MBAR)


Output Files
------------

Results are written under ``{output_dir}/{resname}/``:

``results.txt``
   Human-readable report with LogP and ΔG values per replicate and averaged.

``partition_avg.npy``
   NumPy array of shape ``(n_solvents, 4)`` containing
   ``[dG, dG_err, logP, logP_err]`` for each solvent (averaged over replicates).

``{rep}/{solvent}/{state}/``
   Per-state GROMACS simulation files:

   - ``fep.xvg`` — dH/dλ or u_nk output used for TI/MBAR analysis
   - ``fep.gro`` — final snapshot of the production run
