# %%
import matplotlib
# %matplotlib inline
# print("BACKEND: ", matplotlib.get_backend())
# if matplotlib.get_backend() != "macosx":
#   print("Changing backend to macosx")
#   matplotlib.use('macosx')


# %%
import matplotlib.pyplot as plt
# %matplotlib inline

# matplotlib.use("Agg")

from ase.visualize import view

from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from agox.databases import Database
from agox.environments import Environment
from agox.utils.graph_sorting import Analysis

import pickle
from sklearn.decomposition import PCA


# %%
from agox.models.descriptors.fingerprint import Fingerprint
from agox.models.descriptors import Voronoi
import numpy as np

template = Atoms("", cell=np.eye(3) * 3.567, pbc=True)
diamond = bulk("C", "diamond", a=3.567)  # Lattice constant for diamond cubic carbon
confinement_cell = template.cell.copy()
confinement_corner = np.array([0, 0, 0])
environment = Environment(
    template=template,
    symbols="C8",
    confinement_cell=confinement_cell,
    confinement_corner=confinement_corner,
    box_constraint_pbc=[True, True, True],  # Confinement is periodic in all directions.
)
descriptor = Fingerprint(environment=environment)
graph_descriptor = Voronoi(
    covalent_bond_scale_factor=1.3, n_points=8, angle_from_central_atom=20, environment=environment
)


# %%
# import glob
from chgnet.model import CHGNetCalculator
calc = CHGNetCalculator()
# # for each POSCAR_[0-9] in iteration*/, read in to Atoms() object and attach calculator
# poscar_files = glob.glob("iteration*/POSCAR_[0-9]")
# atoms_list = []
# for poscar_file in poscar_files:
#   atoms = read(poscar_file)
#   atoms.calc =calc
#   atoms_list.append(atoms)


# %%
seed = 0
rlxd_string = "unrlxd"
bin_width = 0.1
target_energy = -9.064090728759766

# %%
for seed in range(1):
  for rlxd_string in ["unrlxd", "rlxd"]:
    structures = read("DOutput/"+rlxd_string+"_structures_seed"+str(seed)+".traj", index=":")
    for structure in structures:
      structure.calc = calc

    # %%
    graph_sorting = Analysis(descriptor=graph_descriptor, directories=["."], sample_size=min(1000, len(structures)))

    # %%
    unique, count = graph_sorting.sort_structures(structures=structures)
    print("Repeats of unique structures: ", count)

    # %%
    # Calculate energies per atom for each unique structure
    energies_per_atom = [structure.get_potential_energy() / len(structure) for structure in unique]
    delta_en_per_atom = np.array(energies_per_atom) - target_energy


    # %%
    if abs( np.min(energies_per_atom) + 9.064090728759766 ) > 5e-2:
      print("Minimum energy per atom is not zero. Check the energy calculation.")
      print("Expected minimum energy per atom: ", -9.064090728759766)
      print("Minimum energy per atom: ", np.min(energies_per_atom))
      print("Difference: ", np.min(energies_per_atom) + 9.064090728759766)
      exit()

    # %%
    # further reduce the delta_en_per_atom and count by the bin width
    delta_en_per_atom_rounded = np.round(delta_en_per_atom / bin_width) * bin_width
    delta_en_per_atom_binned = np.unique(delta_en_per_atom_rounded)
    count_binned = np.zeros_like(delta_en_per_atom_binned)
    for i, de in enumerate(delta_en_per_atom_binned):
      for j, de2 in enumerate(delta_en_per_atom_rounded):
        if abs(de - de2) < 1e-4:
          count_binned[i] += count[j]

    # %%
    # Plot the number of unique structures vs the energies per atom
    plt.figure(figsize=(10, 6))
    plt.bar(delta_en_per_atom_binned, count_binned, width=bin_width)
    plt.xlabel('Energy per Atom (eV/atom)')
    plt.ylabel('Number of structures')
    plt.title('Number of Unique Structures vs Energies per Atom ('+rlxd_string+' seed '+str(seed)+')')
    # plt.show()
    plt.savefig("RAFFLE_unique_structures_vs_energies_per_atom_"+rlxd_string+"_seed"+str(seed)+".png")


    # %%
    pca = PCA(n_components=2)

    # %%
    # save pca model
    # if rlxd_string == "unrlxd" and seed == 0:
    #   pca.fit(np.squeeze([arr for arr in descriptor.get_features(unique)]))
    #   with open("pca_model_"+rlxd_string+"_"+str(seed)+".pkl", "wb") as f:
    #     pickle.dump(pca, f)

    # load pca model
    with open("pca_model_all_unrlxd_0.pkl", "rb") as f:
      pca = pickle.load(f)


    # %%
    X_reduced = pca.transform(np.squeeze([arr for arr in descriptor.get_features(structures)]))
    energies_per_atom_pca = [structure.get_potential_energy() / len(structure) for structure in structures]
    delta_en_per_atom_pca = np.array(energies_per_atom_pca) - target_energy

    # %%
    plt.figure(1, figsize=(8, 6))

    # %%
    plt.figure(figsize=(10, 6))
    plt.scatter(X_reduced[:, 0], X_reduced[:, 1], c=delta_en_per_atom_pca, cmap="viridis")
    plt.xlabel('Principal component 1')
    plt.ylabel('Principal component 2')
    plt.title('PCA of Unique Structures ('+rlxd_string+', seed'+str(seed)+')')
    plt.xlim(-6, 10)
    plt.ylim(-6, 10)
    plt.colorbar(label='Energy above diamond (eV/atom)')
    plt.savefig('RAFFLE_pca_unique_structures_'+rlxd_string+'_seed'+str(seed)+'.png')
