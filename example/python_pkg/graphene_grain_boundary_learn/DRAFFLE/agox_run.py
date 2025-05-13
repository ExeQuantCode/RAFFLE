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
import shephex
from ase import build

@shephex.chain()
def raffle_agox(train=True, method="raffle", index=42):
    # Manually set seed and database-index
    seed = index
    database_index = index
    np.random.seed(seed)

    # Using argparse if e.g. using array-jobs on Slurm to do several independent searches.
    # from argparse import ArgumentParser
    # parser = ArgumentParser()
    # parser.add_argument('-i', '--run_idx', type=int, default=0)
    # args = parser.parse_args()

    # seed = args.run_idx
    # database_index = args.run_idx

    ##############################################################################
    # System & general settings:
    ##############################################################################

    from ase import Atoms
    from ase.io import write
    from ase.calculators.singlepoint import SinglePointCalculator
    from agox.generators.raffle import RaffleGenerator, bounds_to_confinement
    from mace.calculators import mace_mp
    from agox.utils.replica_exchange.priors import get_prior

    chg = get_prior("chg")
    mace = get_prior("mace")

    from ase.build import bulk
    from ase.io import write, read

    template = read("/home/pitfield/RAFFLE/example/python_pkg/graphene_grain_boundary_learn/POSCAR_host_gb")

    graphene = read("/home/pitfield/RAFFLE/example/python_pkg/graphene_grain_boundary_learn/POSCAR_graphene")

    h2 = build.molecule("H2")
    graphene.calc = mace
    C_reference_energy = graphene.get_potential_energy() / len(graphene)
    h2.calc = mace
    H_reference_energy = h2.get_potential_energy() / len(h2)
    element_energies = {
            'C': C_reference_energy,
            'H': H_reference_energy,
        }

    bounds = [[0.20, 0, 0.18], [0.30, 1, 0.33]]
    confinement_corner, confinement_cell  = bounds_to_confinement(bounds, template)

    print(template, confinement_cell, confinement_corner)

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

    if method == "raffle":
        generator = RaffleGenerator(**environment.get_confinement(),
            host=template, 
            database=database,
            environment=environment, 
            element_energies=element_energies,
            order=1, 
            sampler=None,
            train=train, 
            gen_params = {"void": 0.1, "rand": 0.01, "walk": 0.25, "grow": 0.25, "min": 1.0})
    elif method == "rattle":
        if train == "true": exit()
        generator = RattleGenerator(**environment.get_confinement(), sampler=sampler, environment=environment, order=1)

    print("Generator: Setup")
    ##############################################################################
    # Search Settings:
    ##############################################################################

    # Number of steps is very low - should be set higher for a real search!
    relaxer = ParallelRelaxPostprocess(mace, optimizer=BFGS, order=2, optimizer_run_kwargs={"steps": 100, "fmax": 0.05})


    # With the model pre-relax we dont want to take many steps in the real potential!
    # As we are training the model all the data is saved with the store_trajectory argument.
    evaluator = LocalOptimizationEvaluator(
        mace,
        gets={"get_key": "candidates"},
        store_trajectory=False,
        optimizer_run_kwargs={"fmax": 0.05, "steps": 0},
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

    agox = AGOX(generator, minimum_dist, database, evaluator, relaxer, sampler, seed=seed)

    agox.run(N_iterations=1000)

if __name__ == "__main__": 

    walkers = 20
    
    experiments = (raffle_agox(directory='raffle/').
                   zip(train=[True], method=["raffle"]).
                   permute(index=range(25)))
            
    executor = shephex.executor.SlurmExecutor(ntasks=walkers, nodes=1,
            time='0-03:00:00', 
            mem_per_cpu='2G', 
            partition='q48,q40,q36,q28',
            job_name="AGOX-RAFFLE",
            directory="submission", 
            safety_check = False,
            scratch=False)
    
    executor.execute(experiments)