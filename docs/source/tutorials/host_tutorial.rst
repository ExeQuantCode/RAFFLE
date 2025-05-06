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


Preparing the host structure
----------------------------

If a vacuum region has not been set within the host structure (i.e. a region with no atoms between two constituent substructures), the `prepare_host` method can be used to introduce such a region.
The `prepare_host` method will remove atoms from the host structure within a specified distance from the interface location and returns a dictionary of the removed stoichiometry.
The interface location is the location of the interface in the host structure and can be identified using `ARTEMIS <https://artemis-materials.readthedocs.io/en/latest/tutorials/identify_interface_tutorial.html>`_.

Below is an example of how to use both ARTEMIS and RAFFLE to prepare the host structure.

.. code-block:: python

    # Import the required modules
    from artemis.generator import artemis_generator
    artemis_generator = artemis_generator()

    # Return the interface location using ARTEMIS
    location, axis = artemis_generator.get_interface_location( host, return_fractional = True )
    print("Interface found at", location, "along", axis)

    # Prepare the host structure
    missing_stoich = generator.prepare_host(
        interface_location = location,
        depth = 2.0,
        location_as_fractional = True
    )
    print("missing_stoich", missing_stoich)

    # Visualise the prepared host structure
    from ase.visualize import view
    host_modified = generator.get_host()
    view(host_modified)


Defining the bounding box
--------------------------

An optional bounding box can restrict atom placement to a specified region.
The limits are expressed in fractional coordinates relative to the lattice vectors :math:`(\vec{a}, \vec{b}, \vec{c})`.

.. code-block:: python

    # Set the fractional limits of atom position placement
    a_min = 0.0; b_max = 0.0; c_min = 0.3
    a_max = 1.0; b_max = 1.0; c_max = 0.8
    generator.set_bounds( [
        [a_min, b_min, c_min],
        [a_max, b_max, c_max]
    ]  )
