# taster

**taster** is a Python package for computing partition coefficients of Martini 3 small molecules using free energy calculations.

Given a molecule's CG structure and ITP file, taster automates the full workflow:
- Solvation box preparation across multiple solvents
- Parallel FEP/TI simulations via GROMACS
- Free energy estimation using TI or MBAR
- LogP calculation relative to a reference solvent (typically water)
- Convergence and overlap diagnostics (forward/backward convergence, MBAR overlap matrix, TI dhdl plots) for each leg

---

## Requirements

- Python >= 3.10
- GROMACS (2024.3 or greater needed if topologies use angles with restricted bending potentials)
- `numpy >= 1.20`
- `pandas >= 1.0`
- `matplotlib >= 3.5`
- `MDAnalysis >= 2.0`
- `alchemlyb >= 2.0`
- `tqdm >= 4.60`

---

## Installation

```bash
git clone https://github.com/Lp0lp/taster.git
cd taster
pip install -e .
```

---

## Quick Start

### Python API

```python
import taster

# Run the full workflow
df = taster.run_partition_workflow(
    itp='MOL.itp',
    structure='MOL.gro',
    solvents=['water', 'octanol-water_74-26', 'hexadecane', 'chloroform'],
    reference='water',
    reps=3,
    ncores=36,
)

print(df)
```

Or run each step individually:

```python
import taster

# Step 1 — prepare simulation boxes
resname = taster.prepare_partition_setup(
    itp='MOL.itp',
    structure='MOL.gro',
    solvents=['water', 'octanol-water_74-26', 'hexadecane'],
    reps=3,
)

# Step 2 — run TI simulations
taster.run_partitions(resname, solvents=['water', 'octanol-water_74-26', 'hexadecane'],
                      reps=3, ncores=36)

# Step 3 — analyse
df = taster.process_partition(resname, solvents=['octanol-water_74-26', 'hexadecane'],
                               water='water', reps=3)
```

### CLI

```bash
taster --itp MOL.itp --structure MOL.gro \
       --solvents water octanol-water_74-26 hexadecane chloroform \
       --reps 3 --ncores 36 --estimator MBAR
```

See all options:

```bash
taster --help
```

---

## Available Solvents

```python
import taster
print(taster.available_solvents())
```

---

## Output

Results are saved to `./Partitions/<resname>/`:
- `results.txt` — human-readable report of LogP and ΔG values per replicate and average
- `partition_avg.npy` — numpy array of averaged results `(dG, dG_err, logP, logP_err)` per solvent

`process_partition` and `run_partition_workflow` also return a `pandas.DataFrame`:

```
rep  solvent               dG      dG_err   logP    logP_err
1    octanol-water_74-26  -12.3    0.4      2.15    0.07
1    hexadecane           -15.1    0.5      2.64    0.09
...
avg  octanol-water_74-26  -12.3    0.2      2.16    0.04
avg  hexadecane           -15.1    0.3      2.64    0.05
```

### Convergence and Overlap Diagnostics

By default (`diagnostics=True`), `process_partition` also saves convergence/overlap
figures for every replicate and solvent leg (including the reference solvent) under
`./Partitions/<resname>/diagnostics/<rep>/<solvent>/`:
- `mbar_convergence.png` / `ti_convergence.png` — forward/backward convergence of the free energy estimate
- `mbar_overlap_matrix.png` — MBAR overlap matrix between neighbouring lambda states (MBAR only)
- `ti_dhdl.png` — dHdl curve across lambda states (TI only)
- `convergence.log` — convergence analysis log, kept out of the terminal

Pass `diagnostics=False` to skip this and only compute LogP values.

---

## Testing

To run the full test suite from the main folder run:

```bash
pytest
```

This will run all tests in `tests/`. Tests for the analysis module require [alchemtest](https://github.com/alchemistry/alchemtest)

---

## Notes on HPC Usage

By default `ncores` is set to `os.cpu_count()`, which returns the total number of logical CPUs on the machine. On **shared HPC nodes** this will be the node's full core count.

Always pass `ncores` explicitly to match the resources available to your job.

---

## License

LGPLv2.1 License. See [LICENSE](LICENSE) for details.
