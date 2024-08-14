# %%
#Imports MPRester
# from mp_api.client import MPRester
from ase import Atoms
from ase.io import write, read
from pymatgen.io.ase import AseAtomsAdaptor
from mp_api.client import MPRester


# %%
# Personal api key for accessing materials project api
# This is unique to each user
# api_key = ""

# %%
#Gets the ID of materials we want. Ba is all structures with only Ba. Ba-O is all structures with only Ba and O, not only BaO.

mpr = MPRester() # MPRester(api_key)
materials = []

materials.append(mpr.materials.summary.search(chemsys="C", 
                            fields=["material_id","structure", "energy_per_atom", "nsites"]))
materials.append(mpr.materials.summary.search(chemsys="Mg", 
                            fields=["material_id","structure", "energy_per_atom", "nsites"]))
materials.append(mpr.materials.summary.search(chemsys="O", 
                            fields=["material_id","structure", "energy_per_atom", "nsites"]))
materials.append(mpr.materials.summary.search(chemsys="C-Mg", 
                            fields=["material_id","structure", "energy_per_atom", "nsites"]))
materials.append(mpr.materials.summary.search(chemsys="C-O", 
                            fields=["material_id","structure", "energy_per_atom", "nsites"]))
materials.append(mpr.materials.summary.search(chemsys="Mg-O", 
                            fields=["material_id","structure", "energy_per_atom", "nsites"]))


# %%
structures = []
energies = []
nsites = []
for material_set in materials:
    for material in material_set:
        material_id = material.material_id
        structures.append(mpr.get_structure_by_material_id(material_id))
        energies.append(material.energy_per_atom)
        nsites.append(material.nsites)

# %%
all_atoms = []
for structure, energy, nsite in zip(structures, energies, nsites):
    atom = AseAtomsAdaptor.get_atoms(structure)
    atom.info['free_energy'] = energy * nsite
    atom.info['energy'] = energy * nsite
    all_atoms.append(atom)
write("database.xyz", all_atoms, format='extxyz')

# %%



