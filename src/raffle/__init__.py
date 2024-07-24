"""
raffle package

This package provides functionality to interface with a Fortran library,
including a Python wrapper around the Fortran code.
"""

__version__ = '0.2.0'

from raffle.raffle import generator, rw_geom

__all__ = ['__version__', 'generator', 'rw_geom']