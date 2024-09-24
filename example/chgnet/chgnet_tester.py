from ase.io import read, write
# from ase import Atoms
from chgnet.model.dynamics import CHGNetCalculator


# structure = read("test.traj",index=":")
database = read("../example_files/database_perovskites/database.xyz", index=":")

structure = database[0]



calculator = CHGNetCalculator(model=None)

structure.calc = calculator

print(structure.get_potential_energy())

from ase.optimize import BFGS


# trajectory = is the name of the file where the trajectory will be saved
optimizer = BFGS(structure, trajectory = "traje.traj")

optimizer.run(steps=5)


# structure is now the optimized structure
write("structure.vasp", structure)