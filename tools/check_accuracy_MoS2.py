import os
import numpy as np
from pathlib import Path

script_dir = Path(__file__).resolve().parent

## import ASE (Atomic Simulation Environment) modules
from ase import Atoms
from ase.io import read, write


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
database = read("../example/data/MoS2_vasp.xyz", index=":")

# get index of lowest energy structure
lowest_energy_index = np.argmin([
    atoms.get_potential_energy() if atoms.calc is not None else float('inf')
    for atoms in database
])
MoS2_energy_dft = database[lowest_energy_index].get_potential_energy() / len(database[lowest_energy_index])
MoS2_ground_state = database[lowest_energy_index].copy()
MoS2_ground_state.calc = chgnet
MoS2_energy_chgnet = MoS2_ground_state.get_potential_energy() / len(MoS2_ground_state)
MoS2_ground_state.calc = mace
MoS2_energy_mace = MoS2_ground_state.get_potential_energy() / len(MoS2_ground_state)


## Calculate the energies
energies_dft = []
energies_mace = []
energies_chgnet = []
formation_energies_dft = []
formation_energies_mace = []
formation_energies_chgnet = []
for i, atoms in enumerate(database):
    if atoms.calc is None:
        database.remove(atoms)
        continue
    formation_energy_dft = atoms.get_potential_energy() / len(atoms) - MoS2_energy_dft
    formation_energies_dft.append(formation_energy_dft)
    energies_dft.append(atoms.get_potential_energy()/len(atoms))
    atoms.calc = chgnet
    formation_energy_chgnet = atoms.get_potential_energy() / len(atoms) - MoS2_energy_chgnet
    formation_energies_chgnet.append(formation_energy_chgnet)
    energies_chgnet.append(atoms.get_potential_energy()/len(atoms))
    atoms.calc = mace
    formation_energy_mace = atoms.get_potential_energy() / len(atoms) - MoS2_energy_mace
    formation_energies_mace.append(formation_energy_mace)
    energies_mace.append(atoms.get_potential_energy()/len(atoms))


import matplotlib.pyplot as plt

## Write energies to a file
with open("MoS2_energies_comparison.txt", "w") as f:
    f.write("# DFT_Energy_per_atom MACE-MPA-0_Energy_per_atom CHGNet_Energy_per_atom\n")
    for dft_energy, mace_energy, chgnet_energy in zip(energies_dft, energies_mace, energies_chgnet):
        f.write(f"{dft_energy} {mace_energy} {chgnet_energy}\n")

with open("MoS2_formations_comparison.txt", "w") as f:
    f.write("# DFT_Formation_energy_per_atom MACE-MPA-0_Formation_energy_per_atom CHGNet_Formation_energy_per_atom\n")
    for dft_energy, mace_energy, chgnet_energy in zip(formation_energies_dft, formation_energies_mace, formation_energies_chgnet):
        f.write(f"{dft_energy} {mace_energy} {chgnet_energy}\n")

## Plot the energies
plt.figure(figsize=(10, 6))
plt.scatter(energies_dft, energies_mace, c='blue', marker='o', label='MLIP vs DFT')
plt.show()