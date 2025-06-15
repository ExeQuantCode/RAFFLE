.. agox:

=============
AGOX tutorial
=============

RAFFLE is implemented into AGOX as a generator, allowing it to be used in conjunction with other AGOX features such database management and optimisation routines.
The generator is designed to be used with the AGOX framework (:footcite:t:`Christiansen2022AtomisticGlobalOptimization`), which provides a user-friendly way to set up and run RAFFLE-based structure searches.

This tutorial will guide you through the process of setting up and using the RAFFLE generator within AGOX.
The example script can be found in the following directory:

.. code-block:: bash

    raffle/example/python_pkg/agox_runs/Si-Ge_run/agox_run.py

Requirements
------------

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

The RAFFLE generator differs slightly from the standard AGOX generator.
If using `from_host` in RAFFLE, then the host structure must be provided to the generator directly in addition to the environment.
Otherwise, the host structure (without information regarding the system's energy) is provided automatically by the AGOX framework.

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
        element_energies =  {
            'Si': -4.0,
            'Ge': -3.5
        }, # example element energies
        database = database,
        n_structures = 5,
        ...
    )

This sets up the RAFFLE generator to generate 5 structures each iteration, using the host structure and the environment defined earlier.
A more extensive list of arguments specific to the RAFFLE generator can be found at the end of this tutorial in the section :ref:`raffle_generator_arguments`.

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


.. _raffle_generator_arguments:

RAFFLE generator specific arguments
-----------------------------------

The RAFFLE generator has several specific arguments that can be set to control the generation of structures.
The required argument for the RAFFLE generator is:

- ``element_energies``: A dictionary of element energies, taking the form ``{'Si': -4.0, 'Ge': -3.5}``. These are reference energies for the elements in the system, similar to chemical potentials.

More information on this can be found in :ref:`element-energies`.

Other optional arguments that can be set include:

- ``n_structures``: The number of structures to generate in each iteration.
- ``host``: The host structure to use for the generation. This can be an ASE Atoms object.
- ``kBT``: The weighting factor for scaling the importance of different atomic features based on their system's relative energy.
- ``history_len``: The length of the history for tracking change in the generalised descriptor.
- ``width``: The width of the Gaussian functions used in the distribution functions (list of three floats, for 2-, 3-, and 4-body interactions).
- ``sigma``: The standard deviation of the Gaussian functions used in the distribution functions (list of three floats, for 2-, 3-, and 4-body interactions).
- ``cutoff_min``: The minimum cutoff for the Gaussian functions (list of three floats, for 2-, 3-, and 4-body interactions).
- ``cutoff_max``: The maximum cutoff for the Gaussian functions (list of three floats, for 2-, 3-, and 4-body interactions).
- ``radius_distance_tol``: The radius distance tolerance for the element-pair covalent radii (list of four floats, for 3-body min, 3-body max, 4-body min, and 4-body max interactions).
- ``transfer_data``: A list of ASE Atoms objects used to initialise the generalised descriptor.
- ``method_ratio``: The ratio of placement methods to use in the generation (dictionary with keys `void`, `rand`, `walk`, `grow`, and `min`, and values as floats representing the relative importance of each method).
- ``deallocate_systems``: A boolean flag to indicate whether to deallocate the individual distribution functions of systems after they have been combined into the generalised descriptor.
- ``from_host``: A boolean flag to indicate whether to represent the energies with respect to the host structure.
- ``max_walk_attempts``: The maximum number of attempts to place an atom using the `walk` method before giving up.
- ``walk_step_size_coarse``: The initial/coarse step size for the `walk` method.
- ``walk_step_size_fine``: The final/fine step size for the `walk` method.
- ``grid_spacing``: The spacing of the grid used for the `void` and `min` methods.
- ``seed``: A random seed for reproducibility.

Other arguments exist that are not specific to the RAFFLE generator, but are used in the AGOX framework, such as:

- ``database``: The database object to use for storing the generated structures.

Documentation of these arguments can be found in the AGOX documentation: https://agox.gitlab.io/agox/generators/
