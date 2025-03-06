.. quick_guide:

===========
Quick guide
===========


The quick guide is designed to give a brief overview of the steps required to run a RAFFLE calculation.
It is assumed that the user has already installed RAFFLE and has a basic understanding of the command line.
For a more detailed explanation of the parameters and options available, please see the :doc:`Main tutorial page </tutorials/index>` for individual tutorials detailing each step.


RAFFLE is a random structure search package designed primarily for interfaces.
The tutorials will use bulk systems to demonstrate its functionality, but it is recommended to use packages such as `AIRSS <https://airss-docs.github.io>`_ for bulk random structure search.


Single iteration
----------------

Here is a script to run a single iteration of the RAFFLE generator for structure search.

.. code-block:: python

    # Single iteration of RAFFLE structure search
    from ase.io import read, write
    from raffle.generator import raffle_generator

    generator = raffle_generator()

    host = read("host.xyz")
    generator.set_host(host)
    generator.set_grid(grid_spacing=0.1)
    generator.set_bounds([[0, 0, 0.5], [1, 1, 0.75]])
    generator.distributions.set_element_energies(
        { 'C': -9.063733 }
    )

    database = read("database.xyz", index=":")
    generator.distributions.create(database)

    structures, status = generator.generate(
        num_structures = 1,
        stoichiometry = { 'C': 2 },
    )

    write("output.xyz", structures)

This script will generate a single structure with a stoichiometry of two carbon atoms.
The host structure is read from a file, and the grid spacing and bounds are set.
The element energies are set to the energy of carbon, and the database is read from a file.
The generator is then run, and the structures are written to an output file.


Iterative structure search
--------------------------

Here is an example of how to run an iterative structure search with RAFFLE.

.. code-block:: python

    # Iterative RAFFLE structure search
    from ase.io import read, write
    from ase.optimize import FIRE
    from raffle.generator import raffle_generator
    from mace.calculators import mace_mp

    generator = raffle_generator()
    mace = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')

    host = read("host.xyz")
    generator.set_host(host)
    generator.set_grid(grid_spacing=0.1)
    generator.set_bounds([[0, 0, 0.5], [1, 1, 0.75]])
    generator.distributions.set_element_energies(
        { 'C': -9.063733 }
    )

    database = read("database.xyz", index=":")
    generator.distributions.create(database)

    num_structures_old = 0
    for i in range(10):
        generator.generate(
            num_structures = 2,
            stoichiometry = { 'C': 2 }
            calc = calc
        )
        structures, status = generator.get_structures()
        for structure in structures:
            optimiser = FIRE(structure)
            optimiser.run(fmax=0.05)

        generator.update(structures)
        num_structures_old = len(structures)
