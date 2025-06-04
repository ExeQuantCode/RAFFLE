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
def raffle_agox(train=True, n_dist=4, index=42):
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

    from ase.build import bulk
    from ase.io import write
    cell = np.eye(3) * 10 
    template = bulk('C', crystalstructure="diamond", orthorhombic=True)
    template=template.repeat((2,2,4))
    template.calc = mace
    template_copy = template.copy()

    C_reference_energy = template.get_potential_energy() / len(template)

    indices = 0
    
    for i in range(len(template))[::-1]:
        if template[i].position[2] < template.cell[2][2]/2:
            del template[i]
            indices += 1
    
    #for i in range(len(template))[::-1]:
    #    if i > 4:
    #        del template[i]
    #        indices += 1
    #write("template.traj", template)
    print(indices)

    confinement_corner, confinement_cell = bounds_to_confinement([[0, 0, 0.0], [1, 1, 1]],template)

    element_energies = {'C': C_reference_energy}
    species_to_place={"C" : indices}
    
    symbols = ""
    for key in species_to_place.keys(): 
        symbols += key
        symbols += str(species_to_place[key])

    environment = Environment(
        template=template,
        symbols=symbols,
        confinement_cell=confinement_cell,
        confinement_corner=confinement_corner,
    )
    print("Environment: Setup")
    # Database
    db_path = f"../db{database_index}_dist_{n_dist}.db"  # From input argument!
    database = Database(filename=db_path, order=5)
    database.restore_to_memory()

    
    generator = RaffleGenerator(**environment.get_confinement(),
        host=template, 
        database=database,
        environment=environment, 
        species=species_to_place,
        element_energies=element_energies,
        transfer_data=[template_copy],
        order=1, 
        n_dist=n_dist,
        n_structures=1,
        sampler=None,
        train=train,
        seed=seed,
        gen_params = {"void": 0.0, "rand": 0.1, "walk": 0.25, "grow": 0.25, "min": 1.0})
    
    ##############################################################################
    # Search Settings:
    ##############################################################################
    
    # Number of steps is very low - should be set higher for a real search!
    #relaxer = ParallelRelaxPostprocess(mace, 
    #                                   optimizer=BFGS, 
    #                                   order=2, 
    #                                   optimizer_run_kwargs={"steps": 100, "fmax": 0.05},
    #                                   fix_template=True)
    
    
    # With the model pre-relax we dont want to take many steps in the real potential!
    # As we are training the model all the data is saved with the store_trajectory argument.
    evaluator = LocalOptimizationEvaluator(
        mace,
        gets={"get_key": "candidates"},
        store_trajectory=False,
        optimizer_run_kwargs={"fmax": 0.05, "steps": 100},
        order=3,
        number_to_evaluate=1,
        constraints=environment.get_constraints(),
        fix_template = True, 
    )
    from agox.postprocessors.minimum_dist import MinimumDistPostProcess

    c1 = 0.7
    c2 = 3  
    minimum_dist = MinimumDistPostProcess(
        c1=c1,
        c2=c2,
        order=1.5,
    )
    minimum_dist2 = MinimumDistPostProcess(
        c1=0.9,
        gets = {"get_key" : "evaluated_candidates"},
        sets = {"set_key" : "evaluated_candidates"},
        c2=c2,
        order=3.5,
    )
    
    ##############################################################################
    # Let get the show running!
    ##############################################################################

    agox = AGOX(generator, minimum_dist, minimum_dist2, database, evaluator, seed=seed)

    agox.run(N_iterations=1000)

if __name__ == "__main__":

    walkers = 10
    
    experiments = (raffle_agox(directory='raffle/').
                   zip(train=[True]).
                   permute(index=range(10)).
                   permute(n_dist=[2,3,4]))
            
    executor = shephex.executor.SlurmExecutor(ntasks=walkers, nodes=1,
            time='1-00:00:00', 
            mem_per_cpu='10G', 
            partition='q48,q40,q36',
            job_name="AGOX-RAFFLE",
            directory="submission", 
            safety_check = False,
            scratch=False)
    
    executor.execute(experiments)
