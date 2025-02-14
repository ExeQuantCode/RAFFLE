import os
import numpy as np

## import ASE (Atomic Simulation Environment) modules
from ase import Atoms
from ase.io import read, write


## load calculator
calculator = "CHGNet"
match calculator:
    case "CHGNet":
        from chgnet.model.dynamics import CHGNetCalculator
        print("Initialising CHGNet calculator")
        calc = CHGNetCalculator()
        label = "CHGNet"
    case "MACE":
        from mace.calculators import mace_mp
        print("Initialising MACE calculator")
        calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')
        label = "MACE"

## Read the database
print("Reading database")
database = read("../example/data/Si-Ge.xyz", index=":")

Si_energy_dft = -10.644538 / 2 # Si POT
Ge_energy_dft = -8.7301811 / 2 # Ge POT

## Calculate the energies
energies_dft = []
energies_mlp = []
formation_energies_dft = []
formation_energies_mlp = []
for i, atoms in enumerate(database):
    if atoms.calc is None:
        database.remove(atoms)
        continue
    # get number of Si atoms in the structure
    n_Si = atoms.get_atomic_numbers().tolist().count(14)
    n_Ge = atoms.get_atomic_numbers().tolist().count(32)
    if(n_Si == 0 and n_Ge == 0 or n_Si + n_Ge != len(atoms)):
        print(f"Skipping structure {i} as it has no Si or Ge atoms")
        continue
    # if(n_Si == 0 or n_Ge == 0):
    #     print(f"Skipping structure {i} as it has only one type of atom")
    #     continue
    formation_energy_dft = ( atoms.get_potential_energy() - n_Si * Si_energy_dft - n_Ge * Ge_energy_dft ) / len(atoms)
    formation_energies_dft.append(formation_energy_dft)
    energies_dft.append(atoms.get_potential_energy()/len(atoms))
    atoms.calc = mace_calc
    formation_energy_mace = ( atoms.get_potential_energy() - n_Si * Si_energy_dft - n_Ge * Ge_energy_dft ) / len(atoms)
    formation_energies_mlp.append(formation_energy_mace)
    energies_mlp.append(atoms.get_potential_energy()/len(atoms))
    print(n_Si, n_Ge, formation_energy_dft, formation_energy_mace, atoms.get_potential_energy()/len(atoms))
    # if energies_mlp[-1] - energies_dft[-1] > 3e-1 and energies_mlp[-1] < -7.9:
    #     print(f"Energy difference for structure {i} is {energies_mlp[-1] - energies_dft[-1]}, energy_mace: {energies_mlp[-1]}")
    #     view(atoms)


import matplotlib.pyplot as plt

## Write energies to a file
with open("Si-Ge_energies_comparison.txt", "w") as f:
    f.write("# DFT_Energy_per_atom MACE_Energy_per_atom\n")
    for dft_energy, mace_energy in zip(energies_dft, energies_mlp):
        f.write(f"{dft_energy} {mace_energy}\n")

with open("Si-Ge_formations_comparison.txt", "w") as f:
    f.write("# DFT_Formation_energy_per_atom "+label+"_Formation_energy_per_atom\n")
    for dft_energy, mace_energy in zip(formation_energies_dft, formation_energies_mlp):
        f.write(f"{dft_energy} {mace_energy}\n")

## Plotting the energies
plt.figure(figsize=(10, 6))
plt.scatter(energies_dft, energies_mlp, c='blue', marker='o', label=label+' vs DFT')
plt.scatter(formation_energies_dft, formation_energies_mlp, c='blue', marker='o', label=label+' vs DFT')
plt.show()