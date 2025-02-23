# This script demonstrates how to use the raffle generator to generate BaTiO3 perovskite structures

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
host.calc = calculator
print("host energy: ", host.get_potential_energy())

generator.set_host(host)
print("Host read")

# generate bulk BaTiO3 and get its energy
reference_bulk = Atoms("BaTiO3", positions=[[0.0, 0.0, 0.0], 
                                        [2.005, 2.005, 2.005],
                                        [2.005, 2.005, 0.000],
                                        [2.005, 0.000, 2.005],
                                        [0.000, 2.005, 2.005]
                                    ], cell=[4.01, 4.01, 4.01], pbc=True)
reference_bulk.calc = calculator

# generate the bulks of the individual elements and get their energies
Ba_reference = Atoms('Ba4', positions=[
    [1.2746913321062139, 0.0, 1.797511295],
    [3.7591171278937865, 0.0, 5.392533885000001],
    [3.7915955621062145, 3.59177373, 1.7975112950000005],
    [1.242212897893787, 3.59177373, 5.392533885000001]
    ], cell=[5.03380846, 7.18354746, 7.19004518], pbc=True)
Ba_reference.calc = calculator
Ba_reference_energy = Ba_reference.get_potential_energy() / len(Ba_reference)
Ti_reference = Atoms("Ti3", positions=[
    [0.0, 0.0, 0.0],
    [2.2836874152833575, 1.3184875439588069, 1.413122135],
    [2.283687415283358, -1.3184875439588073, 1.413122135]
    ], cell=[[2.2836874152833575, -3.955462631876421, 0.0], [2.2836874152833575, 3.955462631876421, 0.0], [0.0, 0.0, 2.82624427]], pbc=True)
Ti_reference.calc = calculator
Ti_reference_energy = Ti_reference.get_potential_energy() / len(Ti_reference)
O_reference = Atoms("O2", positions=[
    [0.0, 0.0, 0.0],
    [1.23, 0.0, 0.0]], cell=[10, 10, 10], pbc=False)
O_reference.calc = calculator
O_reference_energy = O_reference.get_potential_energy() / len(O_reference)

# print energies
print(f"Ba_reference_energy: {Ba_reference_energy}")
print(f"Ti_reference_energy: {Ti_reference_energy}")
print(f"O_reference_energy: {O_reference_energy}")
print(f"BaTiO3_reference_energy: {reference_bulk.get_potential_energy()}")

# set the element reference energies
print("Setting element energies")
generator.distributions.set_element_energies(
    {
        'Ba': Ba_reference_energy,
        'Ti': Ti_reference_energy,
        'O': O_reference_energy
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

# set the distribution function widths (2-body, 3-body, 4-body)
print("Reading database")
use_database = False
if use_database:
    database = read("database.xyz", index=":")
    for i, atoms in enumerate(database):
        # reset energy to use CHGNet
        atoms.calc = calculator
else:
    database = []
    database.append(reference_bulk)
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
print("Setting bins (discretisation of host cell)")
generator.set_grid(grid_spacing=0.1, grid_offset=[0.0, 0.0, 0.0])
print(generator.grid)

# set the stoichiometry for the structures to be generated
print("Setting stoichiometry to insert")
stoich_dict = { 'Ba': 1, 'Ti': 1, 'O': 3 }

# generate structures
num_structures_old = 0
optimise_structure = False
for iter in range(1):
    print(f"Iteration {iter}")
    print("Generating...")
    # this is the main function to generate structures
    generator.generate(num_structures=1, stoichiometry=stoich_dict, seed=0+iter, verbose=1, method_ratio={"void":0.001, "walk":0.0, "min":1.0})
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
