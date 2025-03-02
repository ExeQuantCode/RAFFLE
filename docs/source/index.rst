======
RAFFLE
======

RAFFLE (pseudoRandom Approach For Finding Local Energetic minima) is a Python and Fortran package for structure prediction applied to interfaces.
RAFFLE can be utilised as a Python package, a Fortran library, or a standalone Fortran executable.
The Python package provides a high-level interface to the Fortran library, which contains the core functionality.

The Python package interfaces seemlessly with `ASE (Atomic Simulation Environment) <https://wiki.fysik.dtu.dk/ase/>`_, allowing for easy reading, writing, and manipulation of atomic structures.
Although the package comes with a built-in atomic structure reader and writer, it is recommended to use ASE due to its greater functionality and wide-reaching support.

The code is provided freely available under the `GNU General Public License v3.0 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.

An example 

.. code-block:: python

    # A simple example of how to use RAFFLE to generate 10 structures of diamond and write them to a single file
    from ase import Atoms
    from ase.io import write
    from ase.calculators.singlepoint import SinglePointCalculator
    from raffle.generator import raffle_generator
    from mace.calculators import mace_mp

    generator = raffle_generator()

    calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')

    host = Atoms('C', positions=[[0, 0, 0]], cell=[10, 10, 10])
    host.calc = calc
    generator.set_host(host)

    generator.distributions.set_element_energies( { 'C': 0.0 } )
    generator.distributions.create(host)

    num_structures_old = 0
    for i in range(10):
        generator.generate(
            num_structures = 2,
            stoichiometry = { 'C': 7 }
        )
        structures = generator.get_structures(calc)
        generator.distributions.update(structures[num_structures_old:])
        num_structures_old = len(structures)

    structures = generator.get_structures(calc)
    for structure in structures:
        structure.calc = SinglePointCalculator(
            structure,
            energy=structure.get_potential_energy(),
            forces=structure.get_forces()
        )
 
    write('structures.traj', structures)

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   about
   install
   tutorials/index
   Python API <modules>

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`