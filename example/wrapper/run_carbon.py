import sys
# caution: path[0] is reserved for script path (or '' in REPL)
# sys.path.insert(1, '../../build/')
import os

import raffle
from ase import Atoms
from ase.io import read, write

from chgnet.model.dynamics import CHGNetCalculator
from ase.optimize import BFGS


print("Initialising CHGNet calculator")
calculator = CHGNetCalculator(model=None)


print("Initialising raffle generator")
generator = raffle.generator.raffle_generator_type()


print("Reading host")
host = read("../example_files/POSCAR_host_carbon")
host_basis = raffle.rw_geom.basis_type(host)
generator.set_host(host_basis)
print("Host read")

print("Setting element energies")
generator.distributions.set_element_energies(
    {
        'C': -9.0266865
    }
)


print("Reading database")
database = read("../example_files/database_carbon/database.xyz", index=":")
num_database = len(database)
database_basis = raffle.rw_geom.basis_type_xnum_array()
database_basis.allocate(num_database)
for i, atoms in enumerate(database):
    database_basis.items[i].fromase(atoms)


print("Database read")

print("Setting database")
generator.distributions.create(database_basis, deallocate_systems=False)
print("Database set")

print("Printing distributions")
generator.distributions.write("distributions.txt")
generator.distributions.write_2body("df2.txt")
generator.distributions.write_3body("df3.txt")
generator.distributions.write_4body("df4.txt")

print("Checking element energies")
print(generator.distributions.get_element_energies())

print("Checking bond radii")
print(generator.distributions.get_bond_radii())

print("Setting bins (discretisation of host cell)")
generator.bins = [12,12,30]

print("Setting stoichiometry to insert")
stoich_list = raffle.generator.stoichiometry_type_xnum_array()
stoich_list.allocate(1)
stoich_list.items[0].element = 'C'
stoich_list.items[0].num = 10

num_structures_old = 0
for iter in range(10):
    print(f"Iteration {iter}")
    print("Generating...")
    generator.generate(num_structures=1, stoichiometry=stoich_list, seed=0, verbose=0, method_probab={"void":1.0, "walk":1.0, "min":1.0})
    print("Generated")

    print("Getting structures")
    # generated_structures = generator.get_structures()
    print("number of structures supposed to be generated: ", generator.num_structures)
    generated_structures = generator.structures
    print("actual number allocated: ",len(generated_structures))
    print("Got structures")

    # check if directory iteration[iter] exists, if not create it
    iterdir = f"iteration{iter}/"
    if not os.path.exists(iterdir):
        os.makedirs(iterdir)

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
        # optimizer = BFGS(atoms, trajectory = "traje.traj")
        # optimizer.run(fmax=0.05)
        atoms.get_potential_energy()
        print(f"Structure {inew} energy: {atoms.get_potential_energy()}")
        print(f"Structure {inew} optimised")
        structures_rlxd.items[inew].fromase(atoms)
        write(iterdir+f"POSCAR_{inew}", atoms)

    print("Updating distributions")
    generator.distributions.update(structures_rlxd, deallocate_systems=False)
    print("Printing distributions")
    generator.distributions.write(iterdir+"distributions.txt")
    generator.distributions.write_2body(iterdir+"df2.txt")
    generator.distributions.write_3body(iterdir+"df3.txt")
    generator.distributions.write_4body(iterdir+"df4.txt")
    num_structures_old = num_structures_new