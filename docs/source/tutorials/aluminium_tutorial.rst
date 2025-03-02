.. aluminium:

==================
Aluminium tutorial
==================

This tutorial will guide you through the process of performing RAFFLE-based structure search for the bulk phases of aluminium.

The tutorial is designed to show how a RAFFLE generator given no prior knowledge can learn bonding and identify known phases.
This is not an expected use-case of RAFFLE due to its application to a bulk system, but still demonstrates the expected workflow and capabilities.

The example script can be found in the following directory:

.. code-block:: bash

    raffle/example/python_pkg/Al_learn/learn.py

We recommend reading through the file and running it to understand the process of learning and generating structures.
However, we will provide a brief overview of the script here.

First, we must import the required packages:

.. code-block:: python

    from ase import Atoms
    from raffle.generator import raffle_generator
    from chgnet.model import CHGNetCalculator
    import numpy as np

Next, we need to set up the RAFFLE generator and the calculator to calculate the energies of the structures.
In this example, we use the CHGNet calculator:

.. code-block:: python

    generator = raffle_generator()

    calc = CHGNetCalculator()

Then, we need to create the host structure.
Here, a set of host cells are generated to represent the multiple potential bulk phases of aluminium.
These are then used to generate the structures.

.. code-block:: python

    crystal_structures = ['orthorhombic', 'hcp']
    hosts = []
    for crystal_structure in crystal_structures:
        for a in np.linspace(3.1, 5.4, num=6):
            atom = build.bulk(
                    name = 'Al',
                    crystalstructure = crystal_structure,
                    a = a, b = a, c = a,
            )
            hosts.append(Atoms(
                'Al',
                positions = [(0, 0, 0)],
                cell = atom.get_cell(),
                pbc = True,
                calculator = calc
            ))

The script then sets parameters for the generator and provides an initial database.
Note, this database is effectively empty, as it only contains a single structure, which is an isolated aluminium atom.

.. code-block:: python

    initial_database = [Atoms('Al', positions=[(0, 0, 0)], cell=[8, 8, 8], pbc=True)]
    initial_database[0].calc = calc
    generator.distributions.create(initial_database)

Finally, the script generates structures using the generator.
The generator is given the host structures.
Finally, the generator is run for each host structure, providing a unique stoichiometry each time and using a custom method ratio.

.. code-block:: python

    for host in hosts:
        generator.set_host(host)
        generator.generate(
            num_structures = 5,
            stoichiometry = { 'Al': num_atoms },
            method_ratio = {"void": 0.5, "rand": 0.001, "walk": 0.5, "grow": 0.0, "min": 1.0},
        )

    structures = generator.get_structures()
    write('structures.traj', structures)
