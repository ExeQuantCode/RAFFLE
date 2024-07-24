# RAFFLE

by Ned Thaddeus Taylor, Joe Pitfield, and Steven Paul Hepplestone

RAFFLE (pseudoRandom Approach For Finding Local Energetic minima) is a package for structural prediction applied to material interfaces. RAFFLE can interface with the [Atomic Simulation Environment (ASE)](https://gitlab.com/ase/ase).

RAFFLE is both a Fortran and a Python library. A standalone Fortran executable is also being developed. The code heavily relies on features of Fortran2008 and above, so there is no backwards compatibility with Fortran95.

## Requirements

- Python 3.11 or later (might work on earlier, have not tested)
- Fortran compiler supporting Fortran 2008 standard or later
- NumPy
- f90wrap
- CMake

Optional:
- ASE

The library bas been developed and tested using the following Fortran compilers:
- gfortran -- gcc 13.2.0

## Installation

To install RAFFLE, the source must be obtained from the git repository. Use the following commands to get started:
```
 git clone https://github.com/nedtaylor/raffle.git
 cd athena
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

For Fortran, CMake is required (fpm in the future, hopefully). The Python library installation can be turned off.
```
mkdir build
cd build
cmake -DBUILD_PYTHON=Off ..
make install
```

This is a quick how-to guide for using RAFFLE.

First, compile using
``` 
make
```

## Using

First, you need to ensure that the following file exists in the directory in which you run RAFFLE:
```
chem.in
```

Each of these files should follow the format found in the current repository. They should each have a header line that starts with "#" and contains the word "element". The example headers should then be followed for filling in data. For the elements.dat, the energy provided can be whatever you want to use as a reference energy. This energy is used for calculating formation energy. The example uses energy/atom of the bulk phase of the element.
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
