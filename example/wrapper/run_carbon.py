# This script demonstrates how to use the raffle generator to generate carbon structures

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
host = read("../example_files/POSCAR_host_graphite_trilayer")
host_basis = raffle.rw_geom.basis_type(host)
write("POSCAR_host", host_basis.toase())
generator.set_host(host_basis)
print("Host read")


# Generate bulk diamond and get its energy
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
# Generate bulk graphite (I think with PBE stacking separation) and get its energy
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

# set the carbon reference energy
is_graphite_reference = False
if is_graphite_reference:
    C_reference_energy = graphite_bulk.get_potential_energy() / 4
    write("POSCAR_ref_bulk", graphite_bulk)
else:
    C_reference_energy = diamond_bulk.get_potential_energy() / 8
    write("POSCAR_ref_bulk", diamond_bulk)
print("Setting element energies")
generator.distributions.set_element_energies(
    {
        'C': C_reference_energy#-9.0266865
    }
)

# set the distribution function widths (2-body, 3-body, 4-body)
generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])

# these are the lower and upper bounds for the bond radii
# the first two values are the lower and upper tolerance on covalent bond radii for 3-body distributions
# the last two values are the lower and upper tolerance on covalent bond radii for 4-body distributions
# the max bondlength cutoff is an upper limit, so you can turn off 3-body and 4-body distributions by 
#    setting all the values above 100.0 (just to be safe)
generator.distributions.set_radius_distance_tol([1.5, 2.5, 3.0, 6.0])

# read in the database of structures to use for generating the distribution functions
print("Reading database")
database = read("../example_files/database_carbon/database.xyz", index=":")
database_basis = raffle.rw_geom.basis_type_xnum_array()

print("Allocating database")
use_database = False
if use_database:
    num_database = len(database)
    database_basis.allocate(num_database)
    for i, atoms in enumerate(database):
        # reset energy to use CHGNet
        atoms.calc = calculator
        print(f"Reading structure {i}")
        database_basis.items[i].fromase(atoms)
else:
    num_database = 1
    database_basis.allocate(num_database)
    # database_basis.items[0].fromase(diamond_bulk)
    database_basis.items[0].fromase(graphite_bulk)
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
print("Getting bins (discretisation of host cell)")
generator.set_grid(grid=[20,20,60], grid_offset=[0.0, 0.0, 0.0])
# generator.set_grid(grid_spacing=0.1, grid_offset=[0.0, 0.0, 0.0])
print(generator.grid)

# set the stoichiometry for the structures to be generated
print("Setting stoichiometry to insert")
stoich_list = raffle.generator.stoichiometry_type_xnum_array()
stoich_list.allocate(1)
stoich_list.items[0].element = 'C'
stoich_list.items[0].num = 35

# generate structures
num_structures_old = 0
optimise_structure = False
for iter in range(1):
    print(f"Iteration {iter}")
    print("Generating...")
    # this is the main function to generate structures
    generator.generate(num_structures=1, stoichiometry=stoich_list, seed=0, verbose=0, method_probab={"void":0.01, "walk":0.0, "min":1.0})
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

    # optimise the structures
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
            optimizer.run(fmax=0.05)
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
