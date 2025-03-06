.. databases:

===================
Building a database
===================

This tutorial will guide you through the process of obtaining a database of structures from the Materials Project database.
The process will use a database of carbon, magnesium, and oxygen structures as an example.
However, this process can be applied to any set of elements in the Materials Project database.

.. note::

    The guide is written for the Materials Project API v0.41.2.
    The guide provided on this page might be out of date.
    If you encounter issues, please refer to the `Materials Project API <https://next-gen.materialsproject.org/api>`_ for the most up-to-date information.
    From our experience, it seems that the API still undergoes regular changes, so it is important to check the website for the most recent information.

A script is provided to obtain a database of carbon structures from the Materials Project database, which can be found in the following directory:

.. code-block:: bash

    raffle/tools/database.py


Requirements
------------

To run the script, you will need to install the following Python packages:

- `pymatgen <https://pymatgen.org/>`_
- `ase <https://wiki.fysik.dtu.dk/ase/>`_
- `mp-api <https://next-gen.materialsproject.org/api>`_

Materials Project API key
-------------------------

Once you have installed the required packages, you need to set up your Materials Project API key.
To do so you must:

1. Create an account on the `Materials Project website <https://next-gen.materialsproject.org/>`_.
2. Obtain your API key from the `API section <https://next-gen.materialsproject.org/api>`_.
3. Run the following command via pymatgen ``pmg config --add PMG_MAPI_KEY <your-api-key>``.

Further details (and up-to-date information) can be found at the following two pages:

- `Materials Project API <https://next-gen.materialsproject.org/api>`_
- `Setting the PMG_MAPI_KEY <https://pymatgen.org/usage.html#setting-the-pmg_mapi_key-in-the-config-file>`_

You might also need to set up a VASP_PSP_DIR environment variable to point to the directory containing the VASP pseudopotentials.
If that is the case, follow the guide provided `here <https://pymatgen.org/installation.html#potcar-setup>`_.

Now we can start the script to obtain the database of structures.
First, we must import the required packages:

.. code-block:: python

    from ase import Atoms
    from ase.io import write
    from pymatgen.io.ase import AseAtomsAdaptor
    from mp_api.client import MPRester

Next, we need to set up the Materials Project API key:

.. code-block:: python

    mpr = MPRester()

Now we can obtain the database of carbon, magnesium, and oxygen structures.

.. code-block:: python

    materials = []
    materials.append(mpr.materials.summary.search(chemsys="C", 
                                fields=["material_id","structure", "energy_per_atom", "nsites"]))
    materials.append(mpr.materials.summary.search(chemsys="Mg", 
                                fields=["material_id","structure", "energy_per_atom", "nsites"]))
    materials.append(mpr.materials.summary.search(chemsys="O", 
                                fields=["material_id","structure", "energy_per_atom", "nsites"]))
    materials.append(mpr.materials.summary.search(chemsys="C-Mg", 
                                fields=["material_id","structure", "energy_per_atom", "nsites"]))
    materials.append(mpr.materials.summary.search(chemsys="C-O", 
                                fields=["material_id","structure", "energy_per_atom", "nsites"]))
    materials.append(mpr.materials.summary.search(chemsys="Mg-O", 
                                fields=["material_id","structure", "energy_per_atom", "nsites"]))

Now that we have the structures, we need to extract the relevant information and save the structures to a file.
This can be done more succinctly with `.traj` files, but for human readability, we will use `.xyz` files (extended XYZ format).

.. code-block:: python

    structures = []
    energies = []
    nsites = []
    for material_set in materials:
        for material in material_set:
            material_id = material.material_id
            structures.append(mpr.get_structure_by_material_id(material_id))
            energies.append(material.energy_per_atom)
            nsites.append(material.nsites)

    all_atoms = []
    for structure, energy, nsite in zip(structures, energies, nsites):
        atom = AseAtomsAdaptor.get_atoms(structure)
        atom.info['free_energy'] = energy * nsite
        atom.info['energy'] = energy * nsite
        all_atoms.append(atom)
    write("database.xyz", all_atoms, format='extxyz')

With this, we now have a database of structures that can be used to initialise the generalised distribution functions in RAFFLE.
The database can be loaded in using the ASE package and then provided as an input to the RAFFLE generator's distributions.

.. code-block:: python

    from ase.io import read
    from raffle.generator import raffle_generator

    atoms = read("database.xyz", index=":")

    generator = raffle_generator()
    generator.distributions.create(atoms)

.. note:: 
    
        You may want to set some of the parameters of the distribution functions.
        If so, this must be done BEFORE calling the `create` method.

We are now ready to generate structures using the database of structures.