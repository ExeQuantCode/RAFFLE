# Python Package Example

This directory contains example Python files demonstrating how to use the raffle package.

## Directory Structure

```
python_pkg/
├── <SYSTEM>          # Example with dataset containing global minimum from which to learn RAFFLE descriptors
    └── run.py        # Script to test placement methods
└── <SYSTEM_learn>    # Example with no prior database, i.e. RAFFLE starts with no prior knowledge
    ├── DRAFFLE       # Directory containing example RAFFLE learning script
        ├── learn.py  # Script to run RSS generating and learning
        └── pca.ipynb # Notebook to plot principal component analysis of RAFFLE results
    └── DRSS          # Directory containing example random structure search script using AGOX
        ├── rss.py    # Script to run RSS
        └── pca.ipynb # Notebook to plot principal component analysis of RSS results
```

## File Descriptions

- `run.py`: Example script highlighting placement method capabilities
- `learn.py`: Example script highlighting learning capabilties of RAFFLE
- `pca.ipynb`: Notebook to generate principal component analysis of generated structures

## Execution Order

To run a RAFFLE example and analyse the results, the following order must be performed.
1. Move to the desired `SYSTEM_learn/DRAFFLE` directory and create a directory to work in:
   ```bash
   cd SYSTEM_learn
   mkdir DOutput
   cd DOutput
   ```

2. Run the `learn.py` script:
   ```bash
   python ../learn.py
   ```

3. Go to parent directory and open the notebook `pca.ipynb`, run all the cells.

NOTE: for the `C_learn/DRSS/` example, the `pca.ipynb` notebook expects that the `C_learn/DRAFFLE/` example script and notebook have been fully run first.
This is because it attempts to use the PCA fit to the RAFFLE results to transform its data.
Doing so enables the two techniques to be properly compared.

The `[un]rlxd_structures_seed0.traj` files are provided (as are the `energies_[un]rlxd+seed0.txt`) as the key outputs from the `learn.py` runs.
For the `rss.py` scripts, the `[un]rlxd_structures_seed0.traj` files are provided.
Whilst other output files are generated during the RAFFLE runs, these are optional files with mostly redundant data.

## Prerequisites

- Python 3.11 or higher
- raffle package installed (`pip install .`)
- `ase` for handling atomic structure data
- `CHGNet`, `MACE-0`, `VASP`, or some other `ase`-compatible calculator
