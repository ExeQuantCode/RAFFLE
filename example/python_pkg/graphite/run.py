# This script demonstrates how to use the raffle generator to generate graphite structures

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
host = read("POSCAR_host_missing_layer")
host.calc = calculator
print("host energy: ", host.get_potential_energy())

generator.set_host(host)
print("Host read")


# generate bulk diamond and get its energy
diamond_bulk = Atoms("C8",
                     positions=[
                         [0.0, 0.0, 1.7803725545451616], 
                         [0.8901862772725809, 0.8901862772725808, 2.6705588318177425],
                         [2.863057429826727e-16, 1.7803725545451616, 1.0901637751067644e-16],
                         [0.8901862772725813, 2.6705588318177425, 0.890186277272581],
                         [1.7803725545451616, 0.0, 1.0901637751067644e-16],
                         [2.6705588318177425, 0.8901862772725808, 0.890186277272581],
                         [1.7803725545451619, 1.7803725545451616, 1.7803725545451619],
                         [2.670558831817743, 2.6705588318177425, 2.670558831817743]
                     ], cell=[
                         3.5607451090903233, 3.5607451090903233, 3.5607451090903233
                     ], pbc=True
)
diamond_bulk.calc = calculator
# generate bulk graphite (I think with PBE stacking separation) and get its energy
graphite_bulk = Atoms("C4", 
                     positions=[
                         [0.0, 0.0, 1.95076825],
                         [0.0, 0.0, 5.85230475],
                         [1.2336456308015413, 0.7122456370278755, 1.95076825],
                         [1.2336456308015415, -0.7122456370278757, 5.85230475]
                     ], cell=[
                         [1.2336456308015413, -2.1367369110836267, 0.0], 
                         [1.2336456308015413,  2.1367369110836267, 0.0],
                         [0.0, 0.0, 7.803073]
                     ], pbc=True
)
graphite_bulk.calc = calculator
defected_graphite_bulk = Atoms("C3", 
                     positions=[
                         [0.0, 0.0, 1.95076825],
                         [0.0, 0.0, 5.85230475],
                         [1.2336456308015413, 0.7122456370278755, 1.95076825]
                     ], cell=[
                         [1.2336456308015413, -2.1367369110836267, 0.0], 
                         [1.2336456308015413,  2.1367369110836267, 0.0],
                         [0.0, 0.0, 7.803073]
                     ], pbc=True
)
defected_graphite_bulk.energy = calculator #graphite_bulk.get_potential_energy()

# set the carbon reference energy
is_graphite_reference = False
if is_graphite_reference:
    C_reference_energy = graphite_bulk.get_potential_energy() / len(graphite_bulk)
    write("POSCAR_ref_bulk", graphite_bulk)
else:
    C_reference_energy = diamond_bulk.get_potential_energy() / len(diamond_bulk)
    write("POSCAR_ref_bulk", diamond_bulk)
print("Setting element energies")
generator.distributions.set_element_energies(
    {
        'C': C_reference_energy
    }
)

# set energy scale
generator.distributions.set_kBT(0.2)
# set the distribution function widths (2-body, 3-body, 4-body)
generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])

# set the radii tolerances for the 3-body and 4-body distributions
#    these are the lower and upper bounds for the bond radii
#    the first two values are the lower and upper tolerance on covalent bond radii for 3-body distributions
#    the last two values are the lower and upper tolerance on covalent bond radii for 4-body distributions
#    the max bondlength cutoff is an upper limit, so you can turn off 3-body and 4-body distributions by 
#       setting all the values above 100.0 (just to be safe)
generator.distributions.set_radius_distance_tol([1.5, 2.5, 3.0, 6.0])

# read in the database of structures to use for generating the distribution functions
print("Reading database")
use_database = False
# By default, the database is not used as it needs to be downloaded manually from the Materials Project.
# A script to download the database can be found in tools/database.py (or a notebook in tools/database.ipynb)
if use_database:
    database = read("database.xyz", index=":")
    for i, atoms in enumerate(database):
        # reset energy to use CHGNet
        atoms.calc = calculator
else:
    database = []
    database.append(graphite_bulk)
    # database.append(defected_graphite_bulk)
    # database.append(diamond_bulk)
print("Database read")

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
print("Getting bins (discretisation of host cell)")
# generator.set_grid(grid=[20,20,60], grid_offset=[0.0, 0.0, 0.0])
generator.set_grid(grid_spacing=0.1, grid_offset=[0.0, 0.0, 0.0])
print(generator.grid)

# set the stoichiometry for the structures to be generated
print("Setting stoichiometry to insert")
stoich_dict = { 'C': 2 }

# generate structures
num_structures_old = 0
optimise_structure = False
for iter in range(1):
    print(f"Iteration {iter}")
    print("Generating...")
    # this is the main function to generate structures
    generator.generate(num_structures=1, stoichiometry=stoich_dict, seed=0+iter, verbose=1, method_ratio={"void":0.0001, "walk":0.0, "min":1.0})
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
