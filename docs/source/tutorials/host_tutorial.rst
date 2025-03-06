.. host:

======================
Setting host structure
======================

This tutorial will detail how to set the host structure for RAFFLE and the optional bounding box.

The host structure is the base structure that will be added to in order to generate the structures.
The bounding box is an optional parameter that can be used to limit the space in which atoms can be placed.


It is recommended to set the host after changing parameters such as the Gaussian parameters and cutoffs (see :doc:`Parameters tutorial </tutorials/parameters_tutorial>`).
However, this is not strictly necessary, as the host can be set at any time.
It is recommended that the host be provided with a calculator, as this will enable calculation of the host structure prior to atom placement.


Follow the parameter tutorial to see how to initialise the generator.

.. code-block:: python

    # Initialise RAFFLE generator
    from raffle.generator import raffle_generator

    generator = raffle_generator()

We shall also initialise the calculator.

.. code-block:: python

    # Set the calculator
    from mace.calculators import mace_mp
    calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')


Defining the host structure
---------------------------

Now we shall initialise an atoms object and set it as the host structure.

.. code-block:: python

    # Set the host structure
    host = Atoms(...)
    host.calc = calc
    generator.set_host(host)


Defining the bounding box
--------------------------

An optional bounding box can restrict atom placement to a specified region.
The limits are expressed in fractional coordinates relative to the lattice vectors :math:`(\vec{a}, \vec{b}, \vec{c})`.

.. code-block:: python

    # Set the fractional limits of atom position placement
    a_min = 0.0; b_max = 0.0; c_min = 0.3
    a_max = 1.0; b_max = 1.0; c_max = 0.8
    generator.set_bounds( [
        [a_min, b_min, c_max],
        [a_max, b_max, c_max]
    ]  )
