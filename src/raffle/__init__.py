"""
raffle package

This package provides functionality to interface with a Fortran library,
including a Python wrapper around the Fortran code.
"""

from importlib.metadata import PackageNotFoundError, version
try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

from .raffle import generator as _generator_class
from .raffle import geom_rw as _geom_rw_class
from .graph_builder import graph_builder as _graph_builder_class
from .nn_fingerprint import NNFingerprint
from .gnn_fingerprint import GNNFingerprint
from .structure_metrics import (
    minimum_image_displacements,
    structure_similarity_rmsd,
    symmetry_aware_displacements,
    wrap_atoms_to_unit_cell,
)
try:
    from .torch_gnn_fingerprint import TorchGNNFingerprint
except Exception:
    TorchGNNFingerprint = None
try:
    from .ensemble_exploration import InferenceEnsemble, InferenceEnsembleConfig
except Exception:
    InferenceEnsemble = None
    InferenceEnsembleConfig = None
# from .raffle import generator


# Use the 'types' module to create simulated 'generator' and 'geom submodules
import types
generator = types.ModuleType('generator')
geom = types.ModuleType('geom')
graph_builder = types.ModuleType('graph_builder')

# Assign the respective class to the simulated 'generator' and 'geom' modules
generator.raffle_generator = _generator_class.raffle_generator
generator.stoichiometry_array = _generator_class.stoichiometry_array

# Assign the class to the simulated 'geom' module
geom.basis_array = _geom_rw_class.basis_array
geom.basis = _geom_rw_class.basis

# Assign the class to the simulated 'graph_builder' module
graph_builder.topology = _graph_builder_class.topology
graph_builder.graph_tensors = _graph_builder_class.graph_tensors


# Add the simulated 'generator' and 'geom' module to the current package
import sys
sys.modules['raffle.generator'] = generator
sys.modules['raffle.geom'] = geom
sys.modules['raffle.graph_builder'] = graph_builder
# Clean up internal imports (remove access to the direct classes)
del _generator_class
del _geom_rw_class
del PackageNotFoundError
del version
del sys
del types
del raffle

__all__ = [
    '__version__',
    'generator',
    'geom',
    'graph_builder',
    'NNFingerprint',
    'GNNFingerprint',
    'TorchGNNFingerprint',
    'minimum_image_displacements',
    'structure_similarity_rmsd',
    'symmetry_aware_displacements',
    'wrap_atoms_to_unit_cell',
]

def __getattr__(name):
    if name == "generator":
        return generator
    elif name == "geom":
        return geom
    elif name == "graph_builder":
        return graph_builder
    raise AttributeError(f"module {__name__} has no attribute {name}")
