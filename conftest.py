import importlib.util
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent


def _load_local_torch_gnn_fingerprint() -> None:
    import raffle as raffle_package

    module_path = REPO_ROOT / "src" / "raffle" / "torch_gnn_fingerprint.py"
    spec = importlib.util.spec_from_file_location("raffle.torch_gnn_fingerprint", module_path)
    if spec is None or spec.loader is None:
        return
    module = importlib.util.module_from_spec(spec)
    sys.modules["raffle.torch_gnn_fingerprint"] = module
    spec.loader.exec_module(module)
    raffle_package.TorchGNNFingerprint = module.TorchGNNFingerprint


_load_local_torch_gnn_fingerprint()

def pytest_addoption(parser):
    parser.addoption(
        "--fortran-compiler",
        action="store",
        default="gfortran",  # Default compiler
        help="Specify the Fortran compiler to use"
    )

def pytest_configure(config):
    # Make the Fortran compiler available globally during tests
    compiler = config.getoption("--fortran-compiler")
    os.environ["FORTRAN_COMPILER"] = compiler
