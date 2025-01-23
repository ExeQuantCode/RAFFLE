# https://agox.gitlab.io/agox/command_line_tools/command_line.html

import matplotlib

matplotlib.use("Agg")


import numpy as np
from ase import Atoms
from ase.io import write

from agox import AGOX
from agox.databases import Database
from agox.environments import Environment
from agox.evaluators import LocalOptimizationEvaluator
from agox.generators import RandomGenerator
from agox.postprocessors import MinimumDistPostProcess


database_index = 0
iteration_directory = "iteration"

##############################################################################
# Calculator
##############################################################################

# from mace.calculators import mace_mp
from chgnet.model import CHGNetCalculator
from ase.build import bulk
import os

# calc = mace_mp()
calc = CHGNetCalculator()

##############################################################################
# System & general settings:
##############################################################################


diamond = bulk("C", "diamond", a=3.567)  # Lattice constant for diamond cubic carbon
diamond.write("diamond.traj")
post1 = MinimumDistPostProcess(order=1.5, c1=0.75)
post2 = MinimumDistPostProcess(order=2.5, c1=0.75)

# Make an empty template with a cubic cell of size 6x6x6 with periodic boundary conditions
# in all directions
a = 3.567
iteration = 0
template = Atoms("", cell=np.eye(3) * a, pbc=True)


# Confinement cell matches the template cell in all dimensions and is placed at
# the origin of the template cell.
confinement_cell = template.cell.copy()
confinement_corner = np.array([0, 0, 0])

environment = Environment(
    template=template,
    symbols="C8",
    confinement_cell=confinement_cell,
    confinement_corner=confinement_corner,
    box_constraint_pbc=[True, True, True],  # Confinement is periodic in all directions.
)

# Database
db_directory = iteration_directory+"{}".format(iteration)+"/"
if not os.path.exists(db_directory):
    os.makedirs(db_directory)
db_path = db_directory+"db{}.db".format(database_index)  # From input argument!

for seed in range(20):
    database = Database(filename=db_path, order=3, initialize=True)

    ##############################################################################
    # Search Settings:
    ##############################################################################

    random_generator = RandomGenerator(**environment.get_confinement(), environment=environment, order=1)

    # Wont relax fully with steps:5 - more realistic setting would be 100+.
    evaluator = LocalOptimizationEvaluator(
        calc,
        gets={"get_key": "candidates"},
        optimizer_run_kwargs={"fmax": 0.05, "steps": 100},
        store_trajectory=True,
        order=2,
        constraints=environment.get_constraints(),
    )

    ##############################################################################
    # Let get the show running!
    ##############################################################################

    agox = AGOX(random_generator, database, evaluator, seed=seed)

    agox.run(N_iterations=1000)

    structure_data = database.get_all_structures_data()
    # check iteration number has changed
    iteration_check = -1
    unrlxd_structures = []
    rlxd_structures = []
    for idx, structure in enumerate(structure_data):
        if structure["iteration"] != iteration_check:
            iteration_check = structure["iteration"]
            unrlxd_structures.append(database.db_to_atoms(structure))
            if idx != 0:
                rlxd_structures.append(database.db_to_atoms(structure_data[idx-1]))
    if len(rlxd_structures) != len(unrlxd_structures):
        rlxd_structures.append(database.db_to_atoms(structure_data[-1]))

    print("Unrelaxed structures", len(unrlxd_structures))
    print("Relaxed structures", len(rlxd_structures))
    write(f"unrlxd_structures_seed{seed}.traj", unrlxd_structures)
    write(f"rlxd_structures_seed{seed}.traj", rlxd_structures)

# database.restore_to_memory()
# structures = database.get_all_candidates() # this is an Atoms object with some more stuff
