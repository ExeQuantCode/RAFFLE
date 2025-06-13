import matplotlib

matplotlib.use("Agg")

from agox import AGOX
from agox.databases import Database
from agox.environments import Environment
from agox.evaluators import LocalOptimizationEvaluator
import shephex
from ase import build

@shephex.chain()
def raffle_agox(index=42, void=0.1, rand=0.01, walk=0.25, grow=0.25, min=1.0):


    # Using argparse if e.g. using array-jobs on Slurm to do several independent searches.
    # from argparse import ArgumentParser
    # parser = ArgumentParser()
    # parser.add_argument('-i', '--run_idx', type=int, default=0)
    # args = parser.parse_args()

    ##############################################################################
    # System & generator settings:
    ##############################################################################

    import numpy as np
    from agox.generators.raffle import RaffleGenerator, bounds_to_confinement
    # from mace.calculators import mace_mp
    from agox.utils.replica_exchange.priors import get_prior
    from ase.io import read

    ## Set seed and database-index
    seed = index
    database_index = index
    np.random.seed(seed)

    ## Set up the MACE calculator
    mace = get_prior("mace-mpa")
    # calc_params = { 'model':  '/FULL/PATH/TO/MODEL/FILE' }
    # mace = mace_mp(**calc_params)

    ## Set up the method ratio dictionary and generate a unique ID
    method_ratio = {
        "void": void,
        "rand": rand,
        "walk": walk,
        "grow": grow,
        "min": min
    }
    # generate unique id from placement dict
    method_id = "_".join(f"{k}{v:g}" for k, v in sorted(method_ratio.items()))
    print("method ratio:", method_ratio)
    print("method_id:", method_id)

    ## Set up the host structure (CHANGE THIS TO YOUR HOST STRUCTURE LOCATION)
    template = read("${HOME}/DAGOX/Dgraphene_grain_boundary/POSCAR_host_gb")
    template.calc = mace
    template.set_pbc(True)

    graphene = read("${HOME}/DAGOX/Dgraphene_grain_boundary/POSCAR_graphene")

    ## Set up the reference energies for C and H
    h2 = build.molecule("H2")
    graphene.calc = mace
    C_reference_energy = graphene.get_potential_energy() / len(graphene)
    h2.calc = mace
    H_reference_energy = h2.get_potential_energy() / len(h2)
    element_energies = {
        'C': C_reference_energy,
        'H': H_reference_energy,
    }
    print("element reference energies:", element_energies)

    ## Set up the species to place in the host structure
    species_to_place={'C' : 15}
    symbols = ""
    for key in species_to_place.keys():
        symbols += key
        symbols += str(species_to_place[key])

    ## Set up the bounding box for the confinement
    confinement_corner, confinement_cell = bounds_to_confinement([[0.2, 0, 0.18], [0.3, 1, 0.33]],template)

    ## Set up the environment
    environment = Environment(
        template=template,
        symbols=symbols,
        confinement_cell=confinement_cell,
        confinement_corner=confinement_corner,
        fix_template = False,
    )

    # Database
    db_path = f"../db{database_index}_{method_id}.db"  # From input argument!
    database = Database(filename=db_path, order=5)
    database.restore_to_memory()

    ## Set up the RaffleGenerator
    generator = RaffleGenerator(
        **environment.get_confinement(),
        host=template,
        database=database,
        environment=environment,
        species=species_to_place,
        element_energies=element_energies,
        transfer_data=[graphene],
        order=1,
        width=[0.04, np.pi/160.0, np.pi/160.0],
        kBT=0.2,
        n_structures=5,
        from_host=False,
        sampler=None,
        train=True,
        seed=seed,
        method_ratio = method_ratio
    )

    ##############################################################################
    # Search Settings:
    ##############################################################################

    ## Set up the structural optimisation evaluator (defaults to ASE's BFGS optimizer)
    # from ase.optimize import FIRE # Uncomment to use FIRE instead of BFGS
    evaluator = LocalOptimizationEvaluator(
        mace,
        # optimizer = FIRE, # Uncomment to use FIRE instead of BFGS
        gets={"get_key": "candidates"},
        store_trajectory=False,
        optimizer_run_kwargs={"fmax": 0.05, "steps": 200},
        order=3,
        number_to_evaluate=5,
        constraints=environment.get_constraints(),
        fix_template = False,
    )

    ## Set up post-processors to filter the generated structures.
    # This checks the minimum distance between atoms in the generated structures.
    from agox.postprocessors.minimum_dist import MinimumDistPostProcess
    c1 = 0.6
    c2 = 5
    minimum_dist = MinimumDistPostProcess(
        c1=c1,
        c2=c2,
        order=1.5,
    )
    minimum_dist2 = MinimumDistPostProcess(
        c1=c1,
        gets = {"get_key" : "evaluated_candidates"},
        sets = {"set_key" : "evaluated_candidates"},
        c2=c2,
        order=3.5,
    )

    ##############################################################################
    # Let get the show running!
    ##############################################################################

    ## Set up the AGOX search workflow
    agox = AGOX(generator, minimum_dist, minimum_dist2, database, evaluator, seed=seed)

    ## Run the AGOX search for N_iterations
    agox.run(N_iterations=40)



if __name__ == "__main__":

    walkers = 16

    ## Set up the experiments to run
    # This is an example of how to set up a series of experiments with different
    # method ratios and different seeds.
    experiments = (raffle_agox(directory='method_ratio/').
                   zip(void=[1.0, 0.1,  0.0, 0.1]).
                   zip(rand=[0.1, 0.01, 1.0, 0.0]).
                   zip(walk=[0.5, 0.25, 0.0, 0.5]).
                   zip(grow=[0.5, 0.25, 0.0, 0.5]).
                   zip(min =[1.0, 1.0,  0.0, 1.0]).
                   permute(index=range(10))
                   )

    ## Uncomment for Slurm execution
    executor = shephex.executor.SlurmExecutor.from_config(
        "basic.json",
        safety_check=False,
        mem_per_cpu='1024',
        time="1-00:00:00",
        ntasks_per_node=walkers
    )

    # ## Uncomment for local execution
    # executor = shephex.executor.LocalExecutor()

    # Submit the experiments to the executor
    executor.execute(experiments)
