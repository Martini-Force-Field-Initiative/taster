# taster

**taster** is a Python package for computing partition coefficients of small molecules using Martini 3 coarse-grained molecular dynamics and thermodynamic integration (TI) free energy calculations.

Given a molecule's CG structure and ITP file, taster automates the full workflow:
- Solvation box preparation across multiple solvents
- Parallel FEP/TI simulations via GROMACS
- Free energy estimation using TI or MBAR
- LogP calculation relative to a reference solvent (typically water)

---

## Requirements

- Python >= 3.9
- GROMACS (tested with 2023.x)
- `numpy`, `pandas`, `MDAnalysis`, `alchemlyb`

---

## Installation

```bash
git clone https://github.com/yourname/taster.git
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

---

## License

SOOOOOOooooonnnn