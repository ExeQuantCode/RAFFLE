# %%
import os
from pathlib import Path
from ase.io import read, write
from ase.visualize import view
from ase.optimize import FIRE
from ase.calculators.singlepoint import SinglePointCalculator
from mace.calculators import mace_mp
from artemis.generator import artemis_generator

script_dir = Path(__file__).resolve().parent

# %%
abrupt = read('SiGe_abrupt_interface_rescaled.vasp')
abrupt.set_pbc(True)

# %%
generator = artemis_generator()

# %%
loc, axis = generator.get_interface_location(abrupt)

# %%
generator.set_shift_method(
    num_shifts = 20
)
# %%
structures = generator.regenerate(abrupt, verbose=1)

# %%
# make directory for the structures
output_dir = script_dir / "SiGe_shifted_structures"
output_dir.mkdir(parents=True, exist_ok=True)
for i, structure in enumerate(structures):
    write(output_dir / f'SiGe_shifted_{i}.vasp', structure, format='vasp', sort=True, direct=True)

exit()
# %%

# check if mace file exists
if not os.path.exists(script_dir  / ".." / ".." / "mace-mpa-0-medium.model"):
    print("MACE-MPA-0 model file not found. Please download the model from the MACE website.")
    print("https://github.com/ACEsuit/mace-foundations/releases/tag/mace_mpa_0")
    exit(1)

# set up the calculator
calc_params = { 'model':  script_dir/ ".." / ".." / "mace-mpa-0-medium.model" }
calc = mace_mp(**calc_params)

# %%
structures = [ abrupt.copy() ]
structures.extend(generator.get_structures(calculator = calc))

for structure in structures:
    structure.calc = calc
    structure.calc = SinglePointCalculator(
        structure,
        energy=structure.get_potential_energy(),
        forces=structure.get_forces()
    )


write('SiGe_shifted_unrlxd.traj', structures)

view(structures)

# %%

for structure in structures:
    structure.calc = calc
    optimizer = FIRE(structure)
    optimizer.run(fmax=0.05, steps=200)
    structure.calc = SinglePointCalculator(
        structure,
        energy=structure.get_potential_energy(),
        forces=structure.get_forces()
    )

write('SiGe_shifted_rlxd.traj', structures)
