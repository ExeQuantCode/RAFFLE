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
database = read("../example/data/carbon.xyz", index=":")

## Calculate the energies
energies_dft = []
energies_mlp = []
for i, atoms in enumerate(database):
    if atoms.calc is None:
        database.remove(atoms)
        continue
    energies_dft.append(atoms.get_potential_energy()/len(atoms))
    atoms.calc = calc
    energies_mlp.append(atoms.get_potential_energy()/len(atoms))
    # if energies_mlp[-1] - energies_dft[-1] > 3e-1 and energies_mlp[-1] < -7.9:
    #     print(f"Energy difference for structure {i} is {energies_mlp[-1] - energies_dft[-1]}, energy_mace: {energies_mlp[-1]}")
    #     view(atoms)


import matplotlib.pyplot as plt

## Write energies to a file
with open("C_energies_comparison.txt", "w") as f:
    f.write("# DFT_Energy_per_atom "+label+"_Energy_per_atom\n")
    for dft_energy, mace_energy in zip(energies_dft, energies_mlp):
        f.write(f"{dft_energy} {mace_energy}\n")

## Plot the energies
plt.figure(figsize=(10, 6))
plt.scatter(energies_dft, energies_mlp, c='blue', marker='o', label=label+' vs DFT')
plt.show()