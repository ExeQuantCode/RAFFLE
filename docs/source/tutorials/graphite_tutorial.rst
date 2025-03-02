.. graphite:

=================
Graphite tutorial
=================

This tutorial will guide you through the process of reconstructing graphite from a defected graphite cell (one missing layer).
This tutorial follows the same structure as the :doc:`Diamond tutorial </tutorials/diamond_tutorial>`.

The tutorial is designed to show how RAFFLE learns from a database and uses this information to generate structures.
This is not an expected use-case of RAFFLE, but merely a demonstration of its capabilities to rebuild a structure from a defected cell.

The example files can be found in the following directory:

.. code-block:: bash

    raffle/example/python_pkg/graphite

First, we need to establish an initial database of structures.
This will be used to initialise the generalised distribution functions, which will inform the placement of atoms in the generated structures.
Here, we detail two routes: 1) using an ASE object of bulk graphite and 2) using the Materials Project database.

Here, we will simply use the ASE object of bulk graphite.
If you wish to use the Materials Project database, please refer to the :doc:`Databases tutorial </tutorials/databases_tutorial>`.

First, we must import the required packages:

.. code-block:: python

    from ase import Atoms
    from raffle.generator import raffle_generator
    from mace.calculators import mace_mp
    import numpy as np

Next, we need to set up the RAFFLE generator and the calculator to calculate the energies of the structures.
In this example, we use the MACE-MP-0 calculator:

.. code-block:: python

    generator = raffle_generator()

    calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')


Then, we need to create the host structure.
This is the base structure that will be added to in order to generate the structures.

.. code-block:: python

    host = Atoms(
            "C6",
            scaled_positions=[
                [0.000000000, 0.000000000, 0.125000000],
                [0.000000000, 0.000000000, 0.375000000],
                [0.000000000, 0.000000000, 0.875000000],
                [0.333333333, 0.666666667, 0.125000000],
                [0.666666667, 0.333333333, 0.375000000],
                [0.666666667, 0.333333333, 0.875000000],
            ], cell=[
                2.4672912616, 2.4672912616, 15.606146000
            ], pbc=True
        )

    host.calc = calc
    generator.set_host(host)

Next, we need to set any variables we want to use for the distributions.
None should be required other than the element energies, which will be used as reference energies to calculate the formation energies of the structures.
However, for completeness, we will set other variables here.
Because there is only one element in this structure, the reference energy is meaningless, so we set it to ``0.0`` here.
This is due to all structures will have the same stoichiometry and thus use the exact same reference energy, so choice of this value is a global shift.

.. code-block:: python

    generator.distributions.set_element_energies( { 'C': 0.0 } )
    generator.distributions.set_kBT(0.2)
    generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])
    generator.distributions.set_radius_distance_tol([1.5, 2.5, 3.0, 6.0])

Now we need to generate the database of structures.
We will provide bulk diamond as the only database entry here.

.. code-block:: python

    database = []
    database.append(
        Atoms("C4", 
            positions=[
                [0.0, 0.0, 1.95076825],
                [0.0, 0.0, 5.85230475],
                [1.2336456308015413, 0.7122456370278755, 1.95076825],
                [1.2336456308015415, -0.7122456370278757, 5.85230475]
            ], cell=[
                [1.2336456308015413, -2.1367369110836267, 0.0], 
                [1.2336456308015413,  2.1367369110836267, 0.0],
                [0.0, 0.0, 7.803073]
            ], pbc=True
        )
    )

This database will now be used to initialise the generalised distribution functions in RAFFLE.

.. code-block:: python

    generator.distributions.create(database)

Finally, we can set the grid on which atom searches are performed (this grid is applied to the host cell).
By default, the grid is generated using a spacing of 0.1 Å.

.. code-block:: python

    generator.set_grid(grid_spacing=0.1, grid_offset=[0.0, 0.0, 0.0])

We are now ready to generate structures using the database of structures.

.. code-block:: python

    num_structures_old = 0
    structures, status = generator.generate(
        num_structures = 1,
        stoichiometry = { 'C': 2 },
        method_ratio = {"void":0.0001, "min":1.0},
        calc = calc
    )

We should now have a structure of layered graphite.
This structure can be visualised using the ASE package.
But this can also be verified energetically.
The generated structure should have double the energy of bulk graphite, found in ```database[0]```.