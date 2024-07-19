import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '../build/')

import build.raffle as raffle
from ase import Atoms
from ase.io import read

# atoms = Atoms('CC', positions=[[0, 0, 0], [1.2, 0, 0]], pbc=True, cell=[2.4, 2.4, 2.4])

print("Initialising raffle generator")
generator = raffle.generator.raffle_generator_type()


print("Reading host")
host = read("POSCAR_host")
host_basis = raffle.rw_geom.bas_type(host)
generator.set_host(host_basis)
print("Host read")

print("Reading database")
database = read("database.xyz", index=":")
num_database = len(database)
database_basis = raffle.rw_geom.bas_type_xnum_array()
database_basis.allocate(num_database)
for i, atoms in enumerate(database):
    database_basis.items[i].fromase(atoms)
print("Database read")

print("Setting database")
generator.distributions.create(database_basis)
print("Database set")

print("Setting bins (discretisation of host cell)")
generator.bins = [50,50,50]

print("Setting stoichiometry to insert")
stoich_list = raffle.generator.stoichiometry_type_xnum_array()
stoich_list.allocate(2)
stoich_list.items[0].element = 'C'
stoich_list.items[0].num = 1
stoich_list.items[1].element = 'Mg'
stoich_list.items[1].num = 2

print("Generating...")
generator.generate(num_structures=2, stoichiometry=stoich_list)
print("Generated")