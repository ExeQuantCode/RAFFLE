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
database = read("../example/data/aluminium.xyz", index=":")

## Compare energy of isolated Al atom (MLPs should not do well with this test)
Al_reference = Atoms('Al', positions=[(0, 0, 0)], cell=[12,12,12], pbc=False)
Al_reference.calc = calc
Al_reference_energy_mlp = Al_reference.get_potential_energy() / len(Al_reference)

Al_reference_energy_dft = -0.19810165
print("Al_reference_energy_mace: ", Al_reference_energy_mlp)
print("Al_reference_energy_dft: ", Al_reference_energy_dft)

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


import matplotlib.pyplot as plt

## Write energies to a file
with open("energies_comparison.txt", "w") as f:
    f.write("# DFT_Energy_per_atom "+label+"_Energy_per_atom\n")
    for dft_energy, mace_energy in zip(energies_dft, energies_mlp):
        f.write(f"{dft_energy} {mace_energy}\n")

## Plotting the energies
plt.figure(figsize=(10, 6))
plt.scatter(energies_dft, energies_mlp, c='blue', marker='o', label=label+' vs DFT')
plt.show()