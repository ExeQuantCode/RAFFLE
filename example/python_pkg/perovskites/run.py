# This script demonstrates how to use the raffle generator to generate perovskite supperlattice structures

# The script reads a host structure from a POSCAR file, and sets it as the host structure for the generator.

# import standard libraries
import os
import numpy as np

# import raffle library
from raffle.generator import raffle_generator

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
generator = raffle_generator()

# read the host structure from a POSCAR file
print("Reading host")
host = read("POSCAR_host")
generator.set_host(host)
print("Host read")

host.calc = calculator
print("host energy: ", host.get_potential_energy())

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
database = read("database.xyz", index=":")


# create the distribution functions
print("Setting database")
generator.distributions.create(database, deallocate_systems=False)
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
stoich_dict = { 'Ba': 4, 'Sr': 4, 'Ti': 8, 'O': 24 }

# generate structures
num_structures_old = 0
optimise_structure = True
for iter in range(20):
    print(f"Iteration {iter}")
    print("Generating...")
    # this is the main function to generate structures
    generator.generate(num_structures=1, stoichiometry=stoich_dict, seed=0+iter, verbose=0, method_ratio={"void":0.001, "walk":0.0, "min":1.0})
    print("Generated")

    print("Getting structures")
    print("number of structures supposed to be generated: ", generator.num_structures)
    generated_structures = generator.get_structures()
    print("actual number allocated: ",len(generated_structures))
    print("Got structures")

    # check if directory iteration[iter] exists, if not create it
    iterdir = f"iteration{iter}/"
    if not os.path.exists(iterdir):
        os.makedirs(iterdir)

    # get energies using MLPs and optimise the structures
    num_structures_new = len(generated_structures)
    for i, atoms in enumerate(generated_structures):
        if(i < num_structures_old):
            continue
        inew = i - num_structures_old
        atoms.calc = calculator
        if optimise_structure:
            optimizer = BFGS(atoms, trajectory = "traje.traj")
            optimizer.run(fmax=0.05)
            print(f"Structure {inew} optimised")
        atoms.get_potential_energy()
        print(f"Structure {inew} energy: {atoms.get_potential_energy()}")
        write(iterdir+f"POSCAR_{inew}", atoms)

    # update the distribution functions
    print("Updating distributions")
    generator.distributions.update(generated_structures[num_structures_old:], deallocate_systems=False)

    # print the new distribution functions to a file
    print("Printing distributions")
    generator.distributions.write(iterdir+"distributions.txt")
    generator.distributions.write_2body(iterdir+"df2.txt")
    generator.distributions.write_3body(iterdir+"df3.txt")
    generator.distributions.write_4body(iterdir+"df4.txt")

    # update the number of structures generated
    num_structures_old = num_structures_new
