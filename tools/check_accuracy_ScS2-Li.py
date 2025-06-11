import os
import numpy as np
from pathlib import Path

script_dir = Path(__file__).resolve().parent

## import ASE (Atomic Simulation Environment) modules
from ase import Atoms
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator


## load calculator
from chgnet.model.dynamics import CHGNetCalculator
print("Initialising CHGNet calculator")
chgnet = CHGNetCalculator()
from mace.calculators import mace_mp
print("Initialising MACE calculator")
calc_params = { 'model': "../example/python_pkg/mace-mpa-0-medium.model" }
mace = mace_mp(**calc_params)
# calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')


## Read the database
print("Reading database")
database = read("../example/data/ScS2-Li.xyz", index=":")
lowest_energy_index = -1
lowest_energy = float('inf')
# Get the lowest energy structure with just the stoichiometry ScS2
for atoms in database:
    if atoms.calc is None:
        continue
    if atoms.get_chemical_formula(empirical=True) == "ScS2" or "S2Sc":
        n_Sc = atoms.get_atomic_numbers().tolist().count(21)
        energy = atoms.get_potential_energy() / n_Sc
        if energy < lowest_energy:
            lowest_energy = energy
            lowest_energy_index = database.index(atoms)
            print(f"Found ScS2 structure with energy {energy} at index {lowest_energy_index}")
    
    if atoms.get_chemical_formula(empirical=True) == "Li":
        print(f"Found Li structure at index {database.index(atoms)}")

n_Sc = database[lowest_energy_index].get_atomic_numbers().tolist().count(21)
Li = read( script_dir / "../example/python_pkg/ScS2-Li_learn/Li.xyz", index=0 )
Li_energy_dft = Li.get_potential_energy() / len(Li)  # Li_AE_GW POT
ScS2_energy_dft = database[lowest_energy_index].get_potential_energy() / n_Sc
ScS2_ground_state = database[lowest_energy_index].copy()

Li_chgnet = Li.copy()
Li_chgnet.calc = chgnet
Li_chgnet.set_constraint(FixAtoms(mask=[True for atom in Li_chgnet]))
ecf = FrechetCellFilter(Li_chgnet)
# structurally optimized Li structure
opt = BFGS(ecf)
opt.run(fmax=0.02)
Li.calc = chgnet
Li_energy_chgnet = Li_chgnet.get_potential_energy() / len(Li_chgnet)

Li_mace = Li.copy()
Li_mace.calc = mace
Li_mace.set_constraint(FixAtoms(mask=[True for atom in Li_mace]))
ecf = FrechetCellFilter(Li_mace)
# structurally optimized Li structure
opt = BFGS(ecf)
opt.run(fmax=0.02)
Li_energy_mace = Li_mace.get_potential_energy() / len(Li_mace)
ScS2_ground_state.calc = chgnet
ScS2_energy_chgnet = ScS2_ground_state.get_potential_energy() / n_Sc
ScS2_ground_state.calc = mace
ScS2_energy_mace = ScS2_ground_state.get_potential_energy() / n_Sc

## Calculate the energies
formation_energies_dft = []
formation_energies_mace = []
formation_energies_chgnet = []
for i, atoms in enumerate(database):
    if atoms.calc is None:
        database.remove(atoms)
        continue
    if not atoms.get_chemical_formula(empirical=True) == "LiS2Sc":
        continue
    # get number of Li atoms in the structure
    n_Li = atoms.get_atomic_numbers().tolist().count(3)
    # get number of Sc atoms in the structure
    n_Sc = atoms.get_atomic_numbers().tolist().count(21)
    formation_energy_dft = ( atoms.get_potential_energy() - n_Li * Li_energy_dft - n_Sc * ScS2_energy_dft ) / n_Sc / n_Li
    formation_energies_dft.append(formation_energy_dft)
    atoms.calc = chgnet
    formation_energy_chgnet = ( atoms.get_potential_energy() - n_Li * Li_energy_chgnet - n_Sc * ScS2_energy_chgnet ) / n_Sc / n_Li
    formation_energies_chgnet.append(formation_energy_chgnet)
    atoms.calc = mace
    formation_energy_mace = ( atoms.get_potential_energy() - n_Li * Li_energy_mace - n_Sc * ScS2_energy_mace ) / n_Sc / n_Li
    formation_energies_mace.append(formation_energy_mace)


## Write energies to a file
with open("ScS2-Li_formations_comparison.txt", "w") as f:
    f.write("# DFT_Formation_energy_per_atom MACE-MPA-0_Formation_energy_per_atom CHGNet_Formation_energy_per_atom\n")
    for dft_energy, mace_energy, chgnet_energy in zip(formation_energies_dft, formation_energies_mace, formation_energies_chgnet):
        f.write(f"{dft_energy} {mace_energy} {chgnet_energy}\n")
