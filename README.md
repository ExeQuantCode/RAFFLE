<p align="center">
<img src="docs/source/RAFFLE_logo_no_background.png" width="250"/>
</p>

[![MIT workflow](https://img.shields.io/badge/License-GPLv3-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html "View GPLv3 license")
[![Latest Release](https://img.shields.io/github/v/release/ExeQuantCode/RAFFLE?sort=semver)](https://github.com/ExeQuantCode/RAFFLE/releases "View on GitHub")
[![Paper](https://img.shields.io/badge/Paper-Phys_Rev_B-blue.svg)](https://link.aps.org/doi/10.1103/PhysRevLett.132.066201)
[![Documentation Status](https://readthedocs.org/projects/raffle-fortran/badge/?version=latest)](https://raffle-fortran.readthedocs.io/en/latest/?badge=latest)
[![FPM](https://img.shields.io/badge/fpm-0.10.1-purple)](https://github.com/fortran-lang/fpm "View Fortran Package Manager")
[![CMAKE](https://img.shields.io/badge/cmake-3.27.7-red)](https://github.com/Kitware/CMake/releases/tag/v3.27.7 "View cmake")
[![GCC compatibility](https://img.shields.io/badge/gcc-14.1.0-green)](https://gcc.gnu.org/gcc-14/ "View GCC")
[![Coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/nedtaylor/48f14ebb5636b54d3813e4b4494903eb/raw/raffle_coverage_main.json)](https://ExeQuantCode.github.io/RAFFLE/ "View coverage report")
[![ReadTheDocs](https://img.shields.io/readthedocs/raffle-fortran)](https://raffle-fortran.readthedocs.io/en/latest/ "RAFFLE ReadTheDocs")


# RAFFLE

by Ned Thaddeus Taylor, Joe Pitfield, and Steven Paul Hepplestone

RAFFLE (pseudoRandom Approach For Finding Local Energetic minima) is a package for structural prediction applied to material interfaces.
RAFFLE can interface with the [Atomic Simulation Environment (ASE)](https://gitlab.com/ase/ase).

RAFFLE is both a Fortran and a Python library, with the option of a Fortran executable.
The code heavily relies on features of Fortran 2018 and above, so there is no backwards compatibility with Fortran95.

## Documentation

Tutorials and documentation are provided on the [docs](http://raffle-fortran.readthedocs.io/) website.
The methodology is detailed in the [Phys Rev B paper](https://link.aps.org/doi/10.1103/PhysRevLett.132.066201).
The software package will be submitted for publication soon.


## Requirements

- Fortran compiler supporting Fortran 2018 standard or later
- fpm or CMake (fpm works only for Fortran installation)

Python-specific installation:

- Python 3.11 or later (might work on earlier, have not tested)
- NumPy.f2py
- f90wrap
- cython
- scikit-build-core
- meson
- make or ninja
- CMake
- ASE (optional)

The library bas been developed and tested using the following Fortran compilers:
- gfortran -- gcc 13.2.0
- gfortran -- gcc 14.1.0

The library is known to not currently work with the intel Fortran compilers.

## Installation

To install RAFFLE, the source must be obtained from the git repository. Use the following commands to get started:
```
 git clone https://github.com/ExeQuantCode/raffle.git
 cd raffle
```


Depending on what language will be used in, installation will fary from this point.

### Python

For Python, the easiest installation is through pip:
```
pip install .
```

Another option is installing it through cmake, which involves:
```
mkdir build
cd build
cmake ..
make install
```

Then, the path to the install directory (`${HOME}/.local/raffle`) needs to be added to the include path. NOTE: this method requires that the user manually installs the `ase`, `numpy` and `f90wrap` modules for Python.

### Fortran

For Fortran, either fpm or cmake are required.

#### fpm

fpm installation is as follows:

```
fpm build --profile release
```

This will install both the Fortran library and the Fortran application for RAFFLE.
The library can then be called from other fpm-built Fortran programs through normal means (usually referencing the location of RAFFLE in the program's own `fpm.toml` file).
The application can be run using **(NOTE: The application is not a priority for development, so is less likely to work)**:
```
fpm run
```

The library can be tested to ensure compilation was successful **(NOTE: Unit tests currently only provide minimal code coverage)**:
```
fpm test --profile release
```

#### cmake

cmake installation is as follows:
```
mkdir build
cd build
cmake [-DBUILD_PYTHON=Off] ..
make install
```
The optional filed (dentoted with `[...]`) can be used to turn off installation of the Python library.
This will build the library in the build/ directory. All library files will then be found in:
```
${HOME}/.local/raffle
```
Inside this directory, the following files will be generated:
```
include/raffle.mod
lib/libraffle.a
```

To check whether RAFFLE has installed correctly and that the compilation works as expected, the following command can be run:
```
ctest
```
This runs the unit tests (found in the `test` directory) to ensure procedures output as expected.

### MacOS

Issues can arise with the built-in C and Fortran compilers when installing on MacOS, particularly for Python installation.
It is recommended to manually install new versions (such as via [Homebrew](https://brew.sh)) and setting the `CC` and `FC` environment variables (i.e. defining the default compilers for the respective language).
We have found the following to work:

```
brew install gcc
brew install gfortran
export CC=$(brew --prefix gfortran)
export FC=$(brew --prefix gcc)
```

Now follow the instructions for the [Python](#python) build methods.


## Examples

After the library has been installed, a set of example programs can be found in the `example` directory (note, the `test` directory is for unit tests to ensure each procedure in the library produces expected outputs after compilation, they are not really needed to be looked at by potential users).

The `example/executable` example uses a shell script to run the Fortran installed application.
It uses a user-editable input file `param.in` in the same directory to change values of the RAFFLE generator and provided database.

The `example/wrapper` directory contains a set of `run_*.py` files that show how RAFFLE can be called using Python to implement it into existing random structure search workflows.
These examples are the recommended ones to run.
To successfully run them, follow the above installation instructions for Python, then go to the `example/wrapper` directory and run one of the scripts.


<!-- 
Next, after the code is compiled, the way to run it is as follows:
```
<PATH TO RAFFLE DIRECTORY>/bin/raffle -f <PARAMETER_FILE>
```

An example parameter file is provided in the repository. This is `param.in`. filename_host is the host structure filename, found in the execution directory. stoichiometry is the number of each element to be added to the structure. elements is a list of the names of chemical elements. It must be the same size as stoichiometry and its size must match the value provided in the num_species tag (this'll be tidied up later).

The code will search for a directory called "database/" and read any enclosing directories for POSCARs and OUTCARs (will be moving to vasprun.xml once the code works). The structure is obtained from the POSCAR, the energy is obtained from the OUTCAR. These structures are used to seed the gvectors (distribution functions). NOTE TO JOE: Ned has added a cutoff function to the 2-body exactly as found in the original Behler and Parrinello paper.

For testing purposes, the code will output the 2-, 3-, and 4-body gvectors to files 2body.txt, 3body.txt, and 4body.txt, respectively, for plotting purposes. For the most part, they look all right right now. Lots more testing needs to be done to ensure they are stable.

The code will then output any generated structures to increment1/strucXXX/POSCAR. Note, if you rerun, the code will likely break as it won't want to write over existing files.

It seems that the void finder works. I don't think that scan or pseudo-random walk work at all (and neither should they as they look for the old directory space to calculated distribution function contributions, instead of using the new gvectors).
-->


## References

If you use this code, please cite our papers:
```text
@article{Pitfield2024PredictingPhaseStability,
  title = {Predicting Phase Stability at Interfaces},
  author = {Pitfield, J. and Taylor, N. T. and Hepplestone, S. P.},
  journal = {Phys. Rev. Lett.},
  volume = {132},
  issue = {6},
  pages = {066201},
  numpages = {8},
  year = {2024},
  month = {Feb},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.132.066201},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.132.066201}
}
```
