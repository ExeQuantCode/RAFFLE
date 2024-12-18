import os

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