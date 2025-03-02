.. si-ge:

========================
Si|Ge interface tutorial
========================

This tutorial will guide you through the process of performing RAFFLE-based structure search for an Si|Ge interface.

The tutorial is designed to show how a RAFFLE generator given no prior knowledge can learn bonding and statistically available interface configurations.
Statistically available configurations are those that are likely to exist due to energetic values close to the ground state, within thermal and entropy conditions.

The example script can be found in the following directory:

.. code-block:: bash

    raffle/example/python_pkg/Si-Ge_learn/learn.py

We recommend reading through the file and running it to understand the process of learning and generating structures.
However, we will provide a brief overview of the script here.

First, we must import the required packages:

.. code-block:: python

    from ase import Atoms
    from raffle.generator import raffle_generator
    from mace.calculators import mace_mp
    import numpy as np

Next, we need to set up the RAFFLE generator and the calculator to calculate the energies of the structures.
In this example, we use the CHGNet calculator:

.. code-block:: python

    generator = raffle_generator()

    calc = mace_mp(model="mace-mpa-0-medium.model")

Note, choice of calculator is important.
The calculator should be valid within and near the chemical environment of the structures being generated; in this case, Si and Ge bonding.
If this is the case, it can accurately identify local minima and provide a good representation of the energy landscape.
For this example, we use the [MACE-MPA-0 calculator](https://github.com/ACEsuit/mace-mp/releases/tag/mace_mpa_0), which, from preliminary testing, has been found to be suitable for Si and Ge bonding.
To use the MACE-MPA-0 calculator, you will need to download the model file from the link provided and place it in the same directory as the script (and install the MACE package using `pip install mace-torch`, version 0.3.10 or later).
The CHGNet calculator was not found to be suitable for Si and Ge interface bonding.

Then, we need to create the host structure.
Here, an abrupt Si|Ge interface is generated.
This is the base structure that will be added to in order to generate the structures.
A vacuum region of 5.54 Å is set up between the Si and Ge regions, in which the atoms will be placed by RAFFLE.

.. code-block:: python

    Si_bulk = build.bulk("Si", crystalstructure="diamond", a=5.43)
    Si_cubic = build.make_supercell(Si_bulk, [[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
    Ge_bulk = build.bulk("Ge", crystalstructure="diamond", a=5.65)
    Ge_cubic = build.make_supercell(Ge_bulk, [[-1, 1, 1], [1, -1, 1], [1, 1, -1]])

    Si_supercell = build.make_supercell(Si_cubic, [[2, 0, 0], [0, 2, 0], [0, 0, 1]])
    Ge_supercell = build.make_supercell(Ge_cubic, [[2, 0, 0], [0, 2, 0], [0, 0, 1]])

    Si_surface = build.surface(Si_supercell, indices=(0, 0, 1), layers=2)
    Ge_surface = build.surface(Ge_supercell, indices=(0, 0, 1), layers=2)

    host = build.stack(Si_surface, Ge_surface, axis=2, distance= 5.43/2 + 5.65/2)
    cell[2, 2] -= 3.8865
    host.set_cell(cell, scale_atoms=False)


The script then sets parameters for the generator and provides an initial database.
Note, this database only contains prior information of Si-Si and Ge-Ge bonding, and no information about Si-Ge bonding.
This is to demonstrate the ability of RAFFLE to learn from scratch.

.. code-block:: python

    Si_bulk.calc = calc
    Ge_bulk.calc = calc
    generator.distributions.set_element_energies(
        {
            'Si': Si_bulk.get_potential_energy() / len(Si_bulk),
            'Ge': Ge_bulk.get_potential_energy() / len(Ge_bulk),
        }
    )

    # set energy scale
    generator.distributions.set_kBT(0.2)

    # set the distribution function widths (2-body, 3-body, 4-body)
    generator.distributions.set_width([0.04, np.pi/160.0, np.pi/160.0])

    # set the initial database
    initial_database = [Si_bulk, Ge_bulk]
    generator.distributions.create(initial_database)

Finally, the script generates structures using the generator.
The generator is given the host structures.
Finally, the generator is run for each host structure, providing a unique stoichiometry each time and using a custom method ratio.
The 

.. code-block:: python

    generator.set_host(host)
    generator.set_bounds([[0, 0, 0.34], [1, 1, 0.52]])
    for iter in range(40):
        # generate the structures
        structures, exit_code = generator.generate(
            num_structures = 5,
            stoichiometry = { 'Si': 16, 'Ge': 16 },
            seed = iter,
            method_ratio = {"void": 0.1, "rand": 0.01, "walk": 0.25, "grow": 0.25, "min": 1.0},
            verbose = 0,
            calc = calc
        )
        generator.distributions.update(structures)

    structures = generator.get_structures()
    write('structures.traj', structures)
