import os
import numpy as np

from mace.calculators import mace_mp
from chgnet.model import CHGNetCalculator
from ase.calculators.singlepoint import SinglePointCalculator

# import raffle library
from raffle.generator import raffle_generator

# import ASE (Atomic Simulation Environment) modules
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from ase import build

from pathlib import Path

script_dir = Path(__file__).resolve().parent

calc_params = { 'model':  script_dir/ ".." / "mace-mpa-0-medium.model" }
calc = mace_mp(**calc_params)

Si_bulk = build.bulk("Si", crystalstructure="diamond", a=5.43)
Si_bulk.calc = calc
Si_reference_energy = Si_bulk.get_potential_energy() / len(Si_bulk)
# Si_cubic = build.make_supercell(Si_bulk, [[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
Ge_bulk = build.bulk("Ge", crystalstructure="diamond", a=5.65)
Ge_bulk.calc = calc
# Ge_cubic = build.make_supercell(Ge_bulk, [[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
Ge_reference_energy = Ge_bulk.get_potential_energy() / len(Ge_bulk)

# set the parameters for the generator
generator = raffle_generator()
generator.distributions.set_element_energies(
    {
        'Si': Si_reference_energy,
        'Ge': Ge_reference_energy,
    }
)

# set energy scale
generator.distributions.set_kBT(0.2)

# set the distribution function widths (2-body, 3-body, 4-body)
generator.distributions.set_width([0.04, np.pi/160.0, np.pi/160.0])
generator.distributions.set_radius_distance_tol([1.5, 2.5, 3.0, 6.0])


# set the initial database
initial_database = [Si_bulk, Ge_bulk]
#generator.distributions.set_bond_radii( {('Si', 'Si'): 2.5, ('Si', 'Ge'): 2.5} )


generator.distributions.create(initial_database, deallocate_systems=False)

f1, f2, f3 = generator.distributions.generate_fingerprint(Si_bulk, atom_index=0)

generator.distributions.write_gdfs("Si_bulk_gdfs.txt")

#print("hi", f1)
print("Adding fingerprint")
generator.distributions.add_fingerprint(
                            ['Si'], [1],
                            0.0, #Si_bulk.get_potential_energy(),
                            f1, f2, f3,
        )

print("DONE")