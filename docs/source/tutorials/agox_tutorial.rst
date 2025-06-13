.. agox:

=============
AGOX tutorial
=============

RAFFLE is implemented into AGOX as a generator, allowing it to be used in conjunction with other AGOX features such database management and optimisation routines.
The generator is designed to be used with the AGOX framework, which provides a user-friendly way to set up and run RAFFLE-based structure searches.

This tutorial will guide you through the process of setting up and using the RAFFLE generator within AGOX.
The example script can be found in the following directory:

.. code-block:: bash

    raffle/example/python_pkg/agox_runs/Si-Ge_run/agox_run.py

This script utilises the `hex` slurm experiment management package.
The script is designed to run on a slurm cluster, but can also be run locally by removing the slurm-specific code.
`hex` can be installed via pip:

.. code-block:: bash

    pip install hex

Documentation for the `hex` package can be found here: https://gitlab.com/Mads-Peter/shephex/

Documentation for the AGOX framework can be found here: https: https://agox.gitlab.io/agox/

Currently, the RAFFLE generator is only available on the `raffle_generator` branch of AGOX.
To use the RAFFLE generator, you need to install the `raffle_generator` branch of AGOX:

.. code-block:: bash

    git clone -b raffle_generator https://gitlab.com/agox/agox.git
    cd agox
    pip install -e .

The RAFFLE generator differs slightly from the standard AGOX generator, as it requires a host structure to be provided directly to the generator.

An environment must be set up for the search, this defines the host structure, the stoichiometry of atoms to be added, and the bounding box for the search.
The bounding box must first be defined in terms of a lower left and upper right corner, which can be done using the `bounds_to_confinement` function.
The `bounds_to_confinement` function is used to convert the bounding box coordinates into the `confinement_corner` and `confinement_cell` parameters required by the `Environment` class in AGOX.

.. code-block:: python

    from agox.generators.raffle import bounds_to_confinement

    confinement_corner, confinement_cell = bounds_to_confinement(
        [
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0]
        ],
        template
    )

The `bounds_to_confinement` function takes a list of two points defining the lower left and upper right corners of the bounding box, and the host structure.

The environment is defined as follows:

.. code-block:: python

    from agox.environments import Environment

    environment = Environment(
        template = template,
        symbols = symbols,
        confinement_cell = confinement_cell,
        confinement_corner = confinement_corner,
        ...
    )

The `template` is the host structure, an ASE atoms object, which can be created using ASE, ARTEMIS, or read in from a file.
The `symbols` is an alphanumeric string defining the stoichiometry of atoms to be added, e.g. 'Si2Ge2' for a 2:2 ratio of Si to Ge.

A database object must be attached to the search, which is used to store the generated structures.
The database can be created using the `Database` class in AGOX in the following way:

.. code-block:: python

    from agox.databases import Database

    db_path = f"../database.db"
    database = Database(filename=db_path, order=5)
    database.restore_to_memory()

The `order` parameter defines the order in which the database is called in the AGOX framework; more can be read about this in the AGOX documentation.

The RAFFLE generator can then be set up using the `RaffleGenerator` class in AGOX:

.. code-block:: python

    from agox.generators.raffle import RaffleGenerator

    generator = RaffleGenerator(
        **environment.get_confinement(),
        host = template,
        database = database,
        environment = environment,
        species = symbols,
        n_structures = 5,
        ...
    )

This sets up the RAFFLE generator to generate 5 structures each iteration, using the host structure and the environment defined earlier.

Evaluators and structure filters can be set up as usual in AGOX.
For example, to set up an evaluator to perform structural optimisation, and a pre- and post-process filter that removes structures with bondlengths less than a certain value, you can use the following code:

.. code-block:: python

    from agox.evaluators import LocalOptimizationEvaluator
    from agox.postprocessors.minimum_dist import MinimumDistPostProcess

    evaluator = LocalOptimizationEvaluator(
        mace,
        gets = {"get_key": "candidates"},
        store_trajectory = False,
        optimizer_run_kwargs = {"fmax": 0.05, "steps": 200},
        order = 3,
        number_to_evaluate = 5,
        constraints = environment.get_constraints(),
        fix_template = False,
    )

    minimum_dist_pre = MinimumDistPostProcess(
        c1 = 0.6,
        c2 = 5,
        order=1.5,
    )
    minimum_dist_post = MinimumDistPostProcess(
        c1 = 0.6,
        gets = {"get_key" : "evaluated_candidates"},
        sets = {"set_key" : "evaluated_candidates"},
        c2 = 5,
        order = 3.5,
    )

Finally, the AGOX search can be set up and run using the `AGOX` class in AGOX:

.. code-block:: python

    agox = AGOX(generator, minimum_dist_pre, minimum_dist_post, database, evaluator, seed=seed)

    ## Run the AGOX search for N_iterations
    agox.run(N_iterations=40)
