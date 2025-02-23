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


# diamond = bulk("Al", "diamond", a=3.567)  # Lattice constant for diamond cubic carbon
# diamond.write("diamond.traj")
# post1 = MinimumDistPostProcess(order=1.5, c1=0.75)
# post2 = MinimumDistPostProcess(order=2.5, c1=0.75)

# # Make an empty template with a cubic cell of size 6x6x6 with periodic boundary conditions
# # in all directions
# a = 3.567
# iteration = 0
# template = Atoms("", cell=np.eye(3) * a, pbc=True)


crystal_structures = [
    'orthorhombic', 'hcp',
]
hosts = []
for crystal_structure in crystal_structures:
    print(f'Crystal structure: {crystal_structure}')
    for a in np.linspace(3.1, 5.4, num=6):
        b = a
        c = a
        atom = bulk(
                name = 'Al',
                crystalstructure = crystal_structure,
                a = a,
                b = b,
                c = c,
        )
        hosts.append(Atoms('Al', positions=[(0, 0, 0)], cell=atom.get_cell(), pbc=True, calculator=calc))
print("number of hosts: ", len(hosts))

iteration = 0
mass = 26.9815385
density = 1.61 # u/A^3
unrlxd_structures = []
rlxd_structures = []
database_index = -1
for template in hosts:
    database_index += 1
    volume = template.get_volume()

    # calculate the number of atoms in the host to achieve the desired density
    num_atoms = round(density * volume / mass) - 1
    if(num_atoms < 1):
        continue
    print(f"Volume: {volume}")
    print(f"Number of atoms: {num_atoms}")

    # Confinement cell matches the template cell in all dimensions and is placed at
    # the origin of the template cell.
    confinement_cell = template.cell.copy()
    confinement_corner = np.array([0, 0, 0])

    environment = Environment(
        template=template,
        symbols=f"Al{num_atoms-1}",
        confinement_cell=confinement_cell,
        confinement_corner=confinement_corner,
        box_constraint_pbc=[True, True, True],  # Confinement is periodic in all directions.
    )

    # Database
    db_directory = iteration_directory+"{}".format(iteration)+"/"
    if not os.path.exists(db_directory):
        os.makedirs(db_directory)
    db_path = db_directory+"db{}.db".format(database_index)  # From input argument!

    for seed in range(1):
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

        agox.run(N_iterations=50)

        structure_data = database.get_all_structures_data()
        # check iteration number has changed
        iteration_check = -1
        for idx, structure in enumerate(structure_data):
            if structure["iteration"] != iteration_check:
                iteration_check = structure["iteration"]
                unrlxd_structures.append(database.db_to_atoms(structure))
                if idx != 0:
                    rlxd_structures.append(database.db_to_atoms(structure_data[idx-1]))
        if len(rlxd_structures) != len(unrlxd_structures):
            rlxd_structures.append(database.db_to_atoms(structure_data[-1]))


for structure in unrlxd_structures:
    structure.calc = None
for structure in rlxd_structures:
    structure.calc = None
print("Unrelaxed structures", len(unrlxd_structures))
print("Relaxed structures", len(rlxd_structures))
write(f"unrlxd_structures_seed{seed}.traj", unrlxd_structures)
write(f"rlxd_structures_seed{seed}.traj", rlxd_structures)

# database.restore_to_memory()
# structures = database.get_all_candidates() # this is an Atoms object with some more stuff
