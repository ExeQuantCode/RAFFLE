# This script demonstrates how to use the raffle generator to generate perovskite supperlattice structures

# The script reads a host structure from a POSCAR file, and sets it as the host structure for the generator.

# import standard libraries
import sys
import os
import numpy as np

# import raffle library
import raffle

# import ASE (Atomic Simulation Environment) modules
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS

# import CHGNet calculator (for MLP energy evaluation)
from chgnet.model.dynamics import CHGNetCalculator

# load CHGNet calculator
print("Initialising CHGNet calculator")
calculator = CHGNetCalculator(model=None)

# set up an instance of the raffle generator
print("Initialising raffle generator")
generator = raffle.generator.raffle_generator_type()

# read the host structure from a POSCAR file
print("Reading host")
host = read("../example_files/POSCAR_host_perovskites")
host_basis = raffle.rw_geom.basis_type(host)
generator.set_host(host_basis)
print("Host read")

# set the element reference energies
print("Setting element energies")
generator.distributions.set_element_energies(
    {
        'Ba': -1.9030744,
        'Sr': -1.5478236,
        'Ti': -7.842818,
        'O': -4.3707458
    }
)

# set the distribution function widths (2-body, 3-body, 4-body)
print("Reading database")
database = read("../example_files/database_perovskites/database.xyz", index=":")
num_database = len(database)
database_basis = raffle.rw_geom.basis_type_xnum_array()
database_basis.allocate(num_database)
for i, atoms in enumerate(database):
    database_basis.items[i].fromase(atoms)
print("Database read")

# create the distribution functions
print("Setting database")
generator.distributions.create(database_basis, deallocate_systems=False)
print("Database set")

# print the distribution functions to a file
print("Printing distributions")
generator.distributions.write("distributions.txt")
generator.distributions.write_2body("df2.txt")
generator.distributions.write_3body("df3.txt")
generator.distributions.write_4body("df4.txt")

# check the element energies and bond radii
print("Checking element energies")
print(generator.distributions.get_element_energies())
print("Checking bond radii")
print(generator.distributions.get_bond_radii())

# set the grid for the host cell
print("Setting bins (discretisation of host cell)")
generator.set_grid(grid=[8,8,128], grid_offset=[0.0, 0.0, 0.0])
print(generator.grid)

# set the stoichiometry for the structures to be generated
print("Setting stoichiometry to insert")
stoich_list = raffle.generator.stoichiometry_type_xnum_array()
stoich_list.allocate(4)
stoich_list.items[0].element = 'Ba'
stoich_list.items[0].num = 4
stoich_list.items[1].element = 'Sr'
stoich_list.items[1].num = 4
stoich_list.items[2].element = 'Ti'
stoich_list.items[2].num = 8
stoich_list.items[3].element = 'O'
stoich_list.items[3].num = 24

# generate structures
num_structures_old = 0
optimise_structure = True
for iter in range(20):
    print(f"Iteration {iter}")
    print("Generating...")
    # this is the main function to generate structures
    generator.generate(num_structures=1, stoichiometry=stoich_list, seed=0+iter, verbose=0, method_probab={"void":0.001, "walk":0.0, "min":1.0})
    print("Generated")

    print("Getting structures")
    print("number of structures supposed to be generated: ", generator.num_structures)
    generated_structures = generator.structures
    print("actual number allocated: ",len(generated_structures))
    print("Got structures")

    # check if directory iteration[iter] exists, if not create it
    iterdir = f"iteration{iter}/"
    if not os.path.exists(iterdir):
        os.makedirs(iterdir)

    # get energies using MLPs and optimise the structures
    print("Converting to ASE")
    num_structures_new = len(generated_structures)
    structures_rlxd = raffle.rw_geom.basis_type_xnum_array()
    structures_rlxd.allocate(num_structures_new - num_structures_old)
    for i, structure in enumerate(generated_structures):
        if(i < num_structures_old):
            continue
        inew = i - num_structures_old
        print(f"Converting structure {i}")
        atoms = structure.toase()
        atoms.calc = calculator
        if optimise_structure:
            optimizer = BFGS(atoms, trajectory = "traje.traj")
            optimizer.run(fmax=0.5)
            print(f"Structure {inew} optimised")
        atoms.get_potential_energy()
        print(f"Structure {inew} energy: {atoms.get_potential_energy()}")
        structures_rlxd.items[inew].fromase(atoms)
        write(iterdir+f"POSCAR_{inew}", atoms)

    # update the distribution functions
    print("Updating distributions")
    generator.distributions.update(structures_rlxd, deallocate_systems=False)

    # print the new distribution functions to a file
    print("Printing distributions")
    generator.distributions.write(iterdir+"distributions.txt")
    generator.distributions.write_2body(iterdir+"df2.txt")
    generator.distributions.write_3body(iterdir+"df3.txt")
    generator.distributions.write_4body(iterdir+"df4.txt")

    # update the number of structures generated
    num_structures_old = num_structures_new
