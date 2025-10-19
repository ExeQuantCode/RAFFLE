from chgnet.model import CHGNetCalculator
from raffle.generator import raffle_generator
from ase import Atoms
from ase.optimize import FIRE
from ase.io import write, read
import numpy as np
import os
from joblib import Parallel, delayed
from pathlib import Path

# set up the calculator
calc_params = {}
calc = CHGNetCalculator()

# read the host
host = read('buckyball.xyz')
host.calc = calc
host_reference_energy = host.get_potential_energy() / len(host)

# get reference energies
Li_reference = Atoms(
    "Li2",
    positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    cell=[3.42, 3.42, 3.42],
    pbc=[True, True, True],
)
Li_reference.calc = calc
Li_reference_energy = Li_reference.get_potential_energy() / len(Li_reference)

generator = raffle_generator()
generator.distributions.set_element_energies(
    {
        'C': 0.0,
        'Li': Li_reference_energy,
    }
)
# set energy scale
generator.distributions.set_kBT(0.4)

# set the distribution function widths (2-body, 3-body, 4-body)
generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])

# set the initial database
initial_database = [host, Li_reference]
generator.distributions.create(initial_database)

generator.set_host(host)
generator.add_bounds(
    shape_bn = "sphere",
    origin = [0.0, 0.0, 0.0],
    lengths = [2.0]
)

# set the parameters for the generator
seed = 0


# generate the structures
generator.generate(
    num_structures = 5,
    stoichiometry = { 'Li': 1 },
    method_ratio = {"void": 0.0, "rand": 1.0, "walk": 0.0, "grow": 0.0, "min": 0.0},
    verbose = 1,
)

# print the number of structures generated
print("Total number of structures generated: ", generator.num_structures)
generated_structures = generator.get_structures(calc)
num_structures_new = len(generated_structures)

# write the generated structures to xyz files
print(generated_structures)
write("generated_structures_sphere.traj", generated_structures)
write("generated_structures_sphere.xyz", generated_structures)
