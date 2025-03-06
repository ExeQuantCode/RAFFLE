.. generating_tutorial:


=====================
Generating structures
=====================


This tutorial will arguments associated with the generating procedure of a RAFFLE generator.

The generating procedure is the process of adding atoms to the host structure to generate new structures.

.. note::
    The host structure must be set before generating structures. See :doc:`Host tutorial </tutorials/host_tutorial>` for more information.


Input arguments
---------------

The generating procedure has several arguments that can be set to control the generation of structures.

The arguments are as follows:

- ``num_structures``: The number of structures to generate.
- ``stoichiometry``: A dictionary of the stoichiometry of atoms to add.
- ``method_ratio``: A dictionary of the ratio of each method to use.
- ``settings_out_file``: The file to write the settings to.
- ``calc``: The calculator to attach to the generated structures.
- ``seed``: The seed for the random number generator.
- ``verbose``: The verbosity level of the generator.


Example
-------

The following example demonstrates how to generate structures using the RAFFLE generator.

.. code-block:: python

    from raffle import Generator
    from mace.calculators import mace_mp

    # Set the host structure
    generator = Generator()
    generator.set_host('host.xyz')

    # Set the generating arguments
    generator.set_generating(num_structures=10,
                              stoichiometry={'Al': 1, 'O': 2},
                              method_ratio={
                                  'min': 5.0,
                                  'walk': 3.0,
                                  'grow': 2.0
                                  'void': 1.0
                                  'rand': 0.1
                              },
                              settings_out_file='settings.json',
                              calc=mace_mp(),
                              seed=42,
                              verbose=1)

    # Generate the structures
    generator.generate()

The above example generates 10 structures with a stoichiometry of 1 Al and 2 O atoms.
The five methods are used with the specified ratios, with the values being renormalised to sum to 1 after the input.

The settings are written to the file ``settings.json`` for future reference to improve reproducibility.
By default, if no file is specified, the settings are not written to a file.

The calculator is attached to the generated structures to calculate the energies of the structures.
The seed is set to 42 to ensure reproducibility of the random number generator.

The verbosity level is set to 2 to provide detailed information about the generation process.
Default verbosity is 0, which provides no output.


Retrieving generated structures
-------------------------------

The generated structures are returned as a list of ASE Atoms objects.

.. code-block:: python

    structures, status = generator.generate()

The ``structures`` variable contains the list of generated structures for that iteration
The ``status`` variable contains the status of the generation process, which can be used to check for errors.
A successful generation will return a status of 0, while an error will return a non-zero status.

Optionally, all structures generated thus far using the generator can be retrieved using the following command:

.. code-block:: python

    all_structures = generator.get_structures()
