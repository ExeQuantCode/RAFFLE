import sys
# caution: path[0] is reserved for script path (or '' in REPL)
# sys.path.insert(1, '../../build/')

import raffle
import ase
from ase import Atoms
from ase.io import read, write

# atoms = Atoms('CC', positions=[[0, 0, 0], [1.2, 0, 0]], pbc=True, cell=[2.4, 2.4, 2.4])

print("Initialising raffle generator")
generator = raffle.generator.raffle_generator_type()


print("Reading host")
host = read("../example_files/POSCAR_host_perovskites")
host_basis = raffle.rw_geom.bas_type(host)
generator.set_host(host_basis)
print("Host read")

print("Setting element energies")
generator.distributions.set_element_energies(
    {
        'Ba': -1.9030744,
        'Sr': -1.5478236,
        'Ti': -7.842818,
        'O': -4.3707458
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


print("Reading database")
database = ase.io.extxyz.read_extxy("../example_files/database_perovskites/database.xyz", index=":")
num_database = len(database)
database_basis = raffle.rw_geom.bas_type_xnum_array()
database_basis.allocate(num_database)
for i, atoms in enumerate(database):
    database_basis.items[i].fromase(atoms)


print("Database read")

print("Setting database")
generator.distributions.create(database_basis)
print("Database set")

print("Printing distributions")
generator.distributions.write("distributions.txt")

print("Checking element energies")
print(generator.distributions.get_element_energies())

print("Checking bond radii")
print(generator.distributions.get_bond_radii())

print("Setting bins (discretisation of host cell)")
generator.bins = [8,8,128]

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

print("Generating...")
generator.generate(num_structures=10, stoichiometry=stoich_list)
print("Generated")

print("Getting structures")
# generated_structures = generator.get_structures()
print("number of structures supposed to be generated: ", generator.num_structures)
generated_structures = generator.structures
print("actual number allocated: ",len(generated_structures))
print("Got structures")

print("Converting to ASE")
for i, structure in enumerate(generated_structures):
    print(f"Converting structure {i}")
    atoms = structure.toase()
    write(f"POSCAR_{i}", atoms)
