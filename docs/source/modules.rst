raffle
======

.. Submodules
.. ----------

.. raffle.raffle module
.. --------------------

.. .. automodule:: raffle.raffle
..    :members:
..    :undoc-members:
..    :show-inheritance:

Module contents
---------------

RAFFLE is a python package for performing structure prediction at interfaces.
The package provides functionality for generating atomic structures by filling in host structures with additional atoms.
The method involves iteratively generating structures and learning the energetically favourable features of the structures (see https://link.aps.org/doi/10.1103/PhysRevLett.132.066201).
The package is built to accommodate energetic and structure data provided by the Atomic Simulation Environment (ASE) package (https://wiki.fysik.dtu.dk/ase/).

Submodules
----------

.. toctree:: 
   :maxdepth: 2

   raffle.generator
   raffle.geom
   raffle.distributions

.. .. automodule:: raffle
..    :members:
..    :undoc-members:
..    :show-inheritance:
