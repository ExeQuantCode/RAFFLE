import sys
# caution: path[0] is reserved for script path (or '' in REPL)
# sys.path.insert(1, '../../build/')
import os

import raffle
import numpy as np
from ase import Atoms
from ase.io import read, write

from chgnet.model.dynamics import CHGNetCalculator
from ase.optimize import BFGS


print("Initialising CHGNet calculator")
calculator = CHGNetCalculator(model=None)


print("Initialising raffle generator")
generator = raffle.generator.raffle_generator_type()


print("Reading host")
host = read("../example_files/POSCAR_host_BaTiO3")
host_basis = raffle.rw_geom.basis_type(host)
generator.set_host(host_basis)
print("Host read")


# Generate bulk diamond and get its energy
reference_bulk = Atoms("BaTiO3", positions=[[0.0, 0.0, 0.0], 
                                        [2.005, 2.005, 2.005],
                                        [2.005, 2.005, 0.000],
                                        [2.005, 0.000, 2.005],
                                        [0.000, 2.005, 2.005]
                                    ], cell=[4.01, 4.01, 4.01], pbc=True)
reference_bulk.calc = calculator


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

print(f"Ba_reference_energy: {Ba_reference_energy}")
print(f"Ti_reference_energy: {Ti_reference_energy}")
print(f"O_reference_energy: {O_reference_energy}")
print(f"BaTiO3_reference_energy: {reference_bulk.get_potential_energy()}")
write("POSCAR_ref_bulk", reference_bulk)


print("Setting element energies")
generator.distributions.set_element_energies(
    {
        'Ba': Ba_reference_energy,
        'Ti': Ti_reference_energy,
        'O': O_reference_energy
    }
)
generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])
# generator.distributions.set_sigma([0.1, 0.025, 0.1])
generator.distributions.set_radius_distance_tol([1.5, 2.5, 3.0, 6.0])
# 3-body WORKS WITH RADIUS DISTANCE TOL 4.0 FOR UPPER !!!
# general issue either with get_bin (MUST CHECK THOROUGHLY!), or smearing closeness of angles


print("Reading database")
database = read("../example_files/database_carbon/database.xyz", index=":")
database_basis = raffle.rw_geom.basis_type_xnum_array()

# num_database = len(database)
# database_basis.allocate(num_database)
# for i, atoms in enumerate(database):
#     # reset energy to use CHGNet
#     atoms.calc = calculator
#     database_basis.items[i].fromase(atoms)

num_database = 1
database_basis.allocate(num_database)
database_basis.items[0].fromase(reference_bulk)


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
generator.bins = [20,20,40]

print("Setting stoichiometry to insert")
stoich_list = raffle.generator.stoichiometry_type_xnum_array()
stoich_list.allocate(3)
stoich_list.items[0].element = 'Ba'
stoich_list.items[0].num = 1
stoich_list.items[1].element = 'Ti'
stoich_list.items[1].num = 1
stoich_list.items[2].element = 'O'
stoich_list.items[2].num = 3

num_structures_old = 0
for iter in range(1):
    print(f"Iteration {iter}")
    print("Generating...")
    generator.generate(num_structures=1, stoichiometry=stoich_list, seed=0, verbose=0, method_probab={"void":0.001, "walk":0.0, "min":1.0})
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
        # print(f"Structure {inew} optimised")
        atoms.get_potential_energy()
        print(f"Structure {inew} energy: {atoms.get_potential_energy()}")
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
