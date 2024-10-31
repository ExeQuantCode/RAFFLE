.. _about:

=====
About
=====


RAFFLE (pseudoRandom Approach For Finding Local Energetic minima) is a package for structural prediction applied to material interfaces.
RAFFLE can interface with the `Atomic Simulation Environment (ASE) <https://gitlab.com/ase/ase>`_.

RAFFLE is both a Fortran and a Python library, with the option of a Fortran executable.
The code heavily relies on features of Fortran 2018 and above, so there is no backwards compatibility with Fortran95.

The library enables users to fill a host structure with additional atoms, where placement of those atoms is determined by a set of placement methods.
These methods are meant to bias towards energetically favourable configurations, whilst still providing a thorough search of the configuration space.