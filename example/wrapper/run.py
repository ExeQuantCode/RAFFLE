# caution: path[0] is reserved for script path (or '' in REPL)
# sys.path.insert(1, '../../build/')

from raffle.generator import raffle_generator
from ase import Atoms
from ase.io import read, write

# atoms = Atoms('CC', positions=[[0, 0, 0], [1.2, 0, 0]], pbc=True, cell=[2.4, 2.4, 2.4])

print("Initialising raffle generator")
generator = raffle_generator()


print("Reading host")
host = read("../example_files/POSCAR_host")
generator.set_host(host)
print("Host read")

print("Setting element energies")
generator.distributions.set_element_energies(
    {
        'C': -9.0266865,
        'Mg': -1.5478236,
        'O': -4.3707458
    }
)

## BOND SETTER
generator.distributions.set_bond_radii(
    {
        ('C', 'C'): 1.5,
        ('C', 'Mg'): 2.0,
        ('C', 'O'): 1.5,
        ('Mg', 'Mg'): 2.0,
        ('Mg', 'O'): 2.0,
        ('O', 'O'): 1.5
    }
)


print("Reading database")
database = read("../example_files/database/database.xyz", index=":")
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
generator.bins = [12,12,30]

print("Setting stoichiometry to insert")
stoich_dict = { 'C': 8, 'Mg': 8 }

print("Generating...")
generator.generate(num_structures=10, stoichiometry=stoich_dict, seed=0, verbose=1, method_probab={"void":1.0, "walk":1.0, "min":1.0})
print("Generated")

print("Getting structures")
# generated_structures = generator.get_structures()
print("number of structures supposed to be generated: ", generator.num_structures)
generated_structures = generator.get_structures()
print("actual number allocated: ",len(generated_structures))
print("Got structures")

print("Converting to ASE")
for i, atoms in enumerate(generated_structures):
    print(f"Converting structure {i}")
    write(f"POSCAR_{i}", atoms)
