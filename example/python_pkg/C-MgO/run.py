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

print("Initialising raffle generator")
generator = raffle_generator()

print("Reading host")
host = read("hosts/POSCAR_1x5_orthorhombic_11.0A_separation")
host.calc = calculator
print("host energy: ", host.get_potential_energy())

generator.set_host(host)
print("Host read")

print("Setting element energies")
generator.distributions.set_element_energies(
    {
        'C':   -9.0266865, # -1.263,
        'Mg':  -1.5478236, # -0.008,
        'O':   -4.3707458, #-1.888
    }
)

# ## BOND SETTER
# generator.distributions.set_bond_radii(
#     {
#         ('C', 'C'): 1.5,
#         ('C', 'Mg'): 2.0,
#         ('C', 'O'): 1.5,
#         ('Mg', 'Mg'): 2.0,
#         ('Mg', 'O'): 2.0,
#         ('O', 'O'): 1.5
#     }
# )

# set energy scale
generator.distributions.set_kBT(0.5)
# set the distribution function widths (2-body, 3-body, 4-body)
generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])

# set the radii tolerances for the 3-body and 4-body distributions
#    these are the lower and upper bounds for the bond radii
#    the first two values are the lower and upper tolerance on covalent bond radii for 3-body distributions
#    the last two values are the lower and upper tolerance on covalent bond radii for 4-body distributions
#    the max bondlength cutoff is an upper limit, so you can turn off 3-body and 4-body distributions by 
#       setting all the values above 100.0 (just to be safe)
generator.distributions.set_radius_distance_tol([1.5, 2.5, 3.0, 6.0])

print("Reading database")
database = read("database.xyz", index=":")
num_database = len(database)

print("Database read")

print("Setting database")
generator.distributions.create(database, deallocate_systems=False)
print("Database set")

print("Printing distributions")
generator.distributions.write("distributions.txt")
generator.distributions.write_2body("df2.txt")
generator.distributions.write_3body("df3.txt")
generator.distributions.write_4body("df4.txt")
# exit(0)

print("Checking element energies")
print(generator.distributions.get_element_energies())

print("Checking bond radii")
print(generator.distributions.get_bond_radii())

print("Setting bins (discretisation of host cell)")
generator.set_grid(grid_spacing=0.2, grid_offset=[0.0, 0.0, 0.0])
print(generator.grid)

print("Setting stoichiometry to insert")
stoich_dict = { 'Mg': 20, 'O': 20 }
generator.distributions.best_energy = 0.0

# generate structures
num_structures_old = 0
optimise_structure = False
for iter in range(10):
    print(f"Iteration {iter}")
    print("Generating...")
    # this is the main function to generate structures
    generator.generate(num_structures=10, stoichiometry=stoich_dict, seed=0+iter, verbose=0, method_ratio={"void":1.0, "walk":1.0, "min":1.0})
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
