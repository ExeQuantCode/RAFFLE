.. graphene_grain_boundary_tutorial:

================================
Graphene grain boundary tutorial
================================

This tutorial will guide you through the process of performing RAFFLE-based structure search for a graphene grain boundary.

The tutorial follows the graphene grain boundary studied by :footcite:t:`Schusteritsch2014PredictingInterfaceStructures`.
The tutorial is designed to show how a RAFFLE generator can identify multiple statistically accessible grain boundary configurations with just the prior knowledge of pristine graphene.
The tutorial also includes a convergence test to determine when the learning process is likely complete.

The example script can be found in the following directory:

.. code-block:: bash

    raffle/example/python_pkg/graphene_grain_boundary/DRAFFLE/learn.py

It is recommended to read through the file and run it to understand the process of learning and generating structures.
It is also recommended to run the script in a subdirectory to avoid cluttering the main directory with files.
However, a brief overview of the script is provided here.

First, we must import the required packages:

.. code-block:: python

    from ase import build
    from ase.io import read
    from ase.constraints import FixAtoms
    from raffle.generator import raffle_generator
    from mace.calculators import mace_mp
    import numpy as np
    from pathlib import Path

    script_dir = Path(__file__).resolve().parent


Next, we need to set up the RAFFLE generator and the calculator to calculate the energies of the structures.
In this example, we use the MACE calculator:

.. code-block:: python

    generator = raffle_generator()
    generator.distributions.set_history_len(8)

    calc = mace_mp(model="mace-mpa-0-medium.model")

.. note::
  Choice of calculator is important.
  The calculator should be valid within and near the chemical environment of the structures being generated.
  If this is the case, it can accurately identify local minima and provide a good representation of the energy landscape.
  This example uses the [MACE-MPA-0 calculator](https://github.com/ACEsuit/mace-mp/releases/tag/mace_mpa_0), which reliably reproduces the reconstructions found in the literature.
  To use the MACE-MPA-0 calculator, you will need to download the model file from the link provided and place it in the same directory as the script (and install the MACE package using `pip install mace-torch`, version 0.3.10 or later).

The host is two sheets of graphene separated by a vacuum region of roughly 3.4 Å (in the literature, the distance is set to 3.0 Å).
The other ends of the sheets are terminated with hydrogen atoms.
The host is read in from a file, but it can also be generated using the ARTEMIS generator.

  .. code-block:: python

    host = read(script_dir  / ".." / "POSCAR_host_gb")
    generator.set_host(host)

We then need to set up the RAFFLE generator by creating the descriptor.

.. code-block:: python

    graphene = read(script_dir / ".." / "POSCAR_graphene")
    h2 = build.molecule("H2")
    graphene.calc = calc
    C_reference_energy = graphene.get_potential_energy() / len(graphene)
    h2.calc = calc
    H_reference_energy = h2.get_potential_energy() / len(h2)
    generator.distributions.set_element_energies(
        {
            'C': C_reference_energy,
            'H': H_reference_energy,
        }
    )

    initial_database = [graphene]
    generator.distributions.create(initial_database)

Finally, the script generates structures using the generator.


.. code-block:: python

    generator.set_host(host)
    generator.set_bounds([[0.20, 0, 0.18], [0.30, 1, 0.33]])
    converged = False
    for iter in range(40):
        # generate the structures
        structures, exit_code = generator.generate(
            num_structures = 5,
            stoichiometry = { 'C': 15 },
            seed = iter,
            method_ratio = {"void": 0.1, "rand": 0.01, "walk": 0.25, "grow": 0.25, "min": 1.0},
            verbose = 0,
            calc = calc
        )

        for structure in structures:
          c = FixAtoms(indices=[atom.index for atom in structure if atom.symbol == 'H'])
          structure.set_constraint(c)

        generator.distributions.update(structures)
        if generator.distributions.is_converged(threshold=1e-3):
            converged = True
            break

    if not converged:
        print("Convergence not reached after 40 iterations.")
    else:
        print("Convergence reached after {} iterations.".format(iter))

    structures = generator.get_structures()
    write('structures.traj', structures)


.. footbibliography::
