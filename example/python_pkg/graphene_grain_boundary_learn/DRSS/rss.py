import matplotlib

matplotlib.use("Agg")

import numpy as np
from ase import Atoms
from ase.optimize import BFGS

from agox import AGOX
from agox.databases import Database
from agox.environments import Environment
from agox.evaluators import LocalOptimizationEvaluator
from agox.generators import RattleGenerator
from agox.models.descriptors import SOAP
from agox.models.GPR import SparseGPR
from agox.models.GPR.kernels import RBF
from agox.models.GPR.kernels import Constant as C
from agox.postprocessors.ray_relax import ParallelRelaxPostprocess
from agox.samplers import MetropolisSampler
from ase import build
from pathlib import Path

script_dir = Path(__file__).resolve().parent

def rss_agox(index=42):
    # Manually set seed and database-index
    seed = index
    database_index = index
    np.random.seed(seed)

    ##############################################################################
    # System & general settings:
    ##############################################################################

    from ase import Atoms
    from ase.io import write
    from ase.calculators.singlepoint import SinglePointCalculator
    from agox.generators import RandomGenerator
    from agox.generators.raffle import bounds_to_confinement
    from mace.calculators import mace_mp
    import os


    # check if mace file exists
    if not os.path.exists(script_dir  / ".." / ".." / "mace-mpa-0-medium.model"):
        print("MACE-MPA-0 model file not found. Please download the model from the MACE website.")
        print("https://github.com/ACEsuit/mace-foundations/releases/tag/mace_mpa_0")
        exit(1)

    # set up the calculator
    # calc_params = {}
    # calc = CHGNetCalculator()
    calc_params = { 'model':  script_dir/ ".." / ".." / "mace-mpa-0-medium.model" }
    mace = mace_mp(**calc_params)

    from ase.build import bulk
    from ase.io import write, read

    template = read(script_dir / ".." / "POSCAR_host_gb")

    bounds = [[0.20, 0, 0.18], [0.30, 1, 0.33]]
    confinement_corner, confinement_cell  = bounds_to_confinement(bounds, template)

    environment = Environment(
        template=template,
        symbols="C15",
        confinement_cell=confinement_cell,
        confinement_corner=confinement_corner,
    )
    print("Environment: Setup")
    # Database
    db_path = "../db{}.db".format(database_index)  # From input argument!
    database = Database(filename=db_path, order=5)
    sampler = MetropolisSampler(temperature=0.25, order=4)

    generator = RandomGenerator(**environment.get_confinement(), environment=environment, order=1)

    print("Generator: Setup")
    ##############################################################################
    # Search Settings:
    ##############################################################################

    # With the model pre-relax we dont want to take many steps in the real potential!
    # As we are training the model all the data is saved with the store_trajectory argument.
    evaluator = LocalOptimizationEvaluator(
        mace,
        gets={"get_key": "candidates"},
        store_trajectory=False,
        optimizer_run_kwargs={"fmax": 0.05, "steps": 200},
        order=3,
        constraints=environment.get_constraints(),
    )
    from agox.postprocessors.minimum_dist import MinimumDistPostProcess

    c1 = 0.8
    c2 = 3  
    minimum_dist = MinimumDistPostProcess(
        c1=c1,
        c2=c2,
        order=2.5,
    )
    ##############################################################################
    # Let get the show running!
    ##############################################################################

    agox = AGOX(generator, minimum_dist, minimum_dist, database, evaluator, seed=seed)

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


if __name__ == "__main__": 
    
    rss_agox()