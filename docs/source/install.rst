.. _install:

============
Installation
============

For the Python library, the easiest method of installation is to install it directly from pip:

.. code-block:: bash

    pip install raffle

or

.. code-block:: bash

    pip install 'raffle[ase]'

Once this is done, RAFFLE is ready to be used.

Alternatively, to install RAFFLE from source, follow the instructions below.


RAFFLE can be installed in one of three ways; as a Python package, as a Fortran library, or as a standalone Fortran executable.
All versions rely on the core Fortran code, with the Python package and standalone executable wrapping this code in a Python and Fortran interface, respectively.

The code is hosted on `GitHub <https://github.com/ExeQuantCode/raffle>`_.

This can be done by cloning the repository:

.. code-block:: bash

    git clone https://github.com/ExeQuantCode/RAFFLE
    cd RAFFLE

.. note::
    Notes on the different versions of RAFFLE:

    - The Python package is the most user-friendly way to use RAFFLE, as it provides a high-level interface to the Fortran code and allows for easy implementation into existing workflows.

    - The Fortran library provides the same functionality as the Python package, but in Fortran instead, enabling slightly more (but not recommended) access to the underlying code.

    - The standalone executable enables the use of RAFFLE without any Python dependencies and no need to code, just editing of a parameter file.


Global requirements
===================

All installation methods require the following dependency:

- Fortran compiler (gfortran>=13.1, not compatible with intel compilers)

Python
======

Requirements
------------

- python (>=3.11)
- `pip <https://pip.pypa.io/en/stable/>`_
- `f90wrap <https://github.com/jameskermode/f90wrap>`_ (>=0.2.14)
- `numpy <https://numpy.org>`_ (>=1.26)
- `meson <https://mesonbuild.com>`_ (>=1.6)
- `cython <https://cython.org>`_ (>=3.0)
- `sckit-build-core <https://scikit-build-core.readthedocs.io/en/latest/>`_ (>=0.11)
- `cmake <https://cmake.org>`_ (>=3.17)
- `ninja <https://ninja-build.org>`_ (>=1.10) or `GNU Make <https://www.gnu.org/software/make/>`_

Optional
--------

- `ASE <https://wiki.fysik.dtu.dk/ase/>`_ (>=3.23)

Installation using pip
-----------------------

The easiest way to install RAFFLE is via pip.
The package is directly available via PyPI, so can be installed without downloading the repository. To do so, run:

.. code-block:: bash

    pip install 'raffle[ase]'

This will install the RAFFLE package and all its dependencies in the default location.
This is the recommended method of installation, as it is the easiest and most straightforward way to get started with RAFFLE.
By default, this will install the parallel version of RAFFLE.
To install the serial version of RAFFLE, include the ``--config-settings`` flag in the pip install command:

.. code-block:: bash

    pip install --upgrade raffle --config-settings="cmake.define.CMAKE_BUILD_TYPE=Serial"


Another option is to install RAFFLE from the source code, which is recommended if you want to use the latest version of RAFFLE or if you want to contribute to the development of RAFFLE.
To do this, you will need to clone the repository from GitHub.

Once the library is cloned, navigate to the root directory of the repository and run:

.. code-block:: bash

    pip install --upgrade .[ase]


The optional ``[ase]`` flag is used to install the ASE interface, which is not required for the core functionality of RAFFLE, but is recommended (and examples rely on it).

To install the serial version of RAFFLE, include the ``--config-settings`` flag in the pip install command:

.. code-block:: bash

    pip install --upgrade . --config-settings="cmake.define.CMAKE_BUILD_TYPE=Serial"

.. note::
    If you are installing RAFFLE on a computer where you do not have root access, you may need to add the ``--user`` flag to the above command.

Depending on your setup, this will install the Python package and all its dependencies in different places.
To find where this has been installed, you can run:

.. code-block:: bash

    pip show raffle

This will show you the location of the installed package, in addition to other information about the package.

Installation using cmake
------------------------

If you would like to install RAFFLE using cmake, you can do so by running the following commands:

.. code-block:: bash

    mkdir build
    cd build
    cmake [-DBUILD_PYTHON=On] -DBUILD_EXECUTABLE=Off ..
    make
    make install

This will build the Python package and install it in the default location.
For Unix systems, this will typically be in:

.. code-block:: bash

    ~/.local/raffle

Like with the pip installation, this will install the parallel version of RAFFLE by default.
To install the serial version of RAFFLE, run the following commands:

.. code-block:: bash

    mkdir build
    cd build
    cmake [-DBUILD_PYTHON=On] -DBUILD_EXECUTABLE=Off -DCMAKE_BUILD_TYPE=Serial ..
    make
    make install


Fortran
=======

Requirements
------------

- `cmake <https://cmake.org>`_ (>=3.17) or `fpm <https://fpm.fortran-lang.org>`_ (>=0.9.0)
- `GNU Make <https://www.gnu.org/software/make/>`_ (if using cmake)


As mentioned, the Fortran library provides the same functionality as the Python package, but in Fortran instead.

To install the Fortran library or executable, the recommended method is to use the Fortran package manager (fpm).
Cmake is also supported.

Installation using fpm
----------------------

To install the Fortran library and the executable using fpm, navigate to the root directory of the repository and run:

.. code-block:: bash

    fpm build
    fpm install

This can also be set up as a dependency in your own fpm project by adding the following to your ``fpm.toml`` file:

.. code-block:: toml

    [dependencies]
    raffle = { git = "https://github.com/ExeQuantCode/RAFFLE" }

The Fortran libray has openmp dependencies.

Installation using cmake
------------------------

To install the Fortran library using cmake, navigate to the root directory of the repository and run:

.. code-block:: bash

    mkdir build
    cd build
    cmake -DBUILD_PYTHON=Off -DBUILD_EXECUTABLE=Off ..
    make
    make install

This will build the Fortran library and install it in the default location (``~/.local/raffle``).

To install the standalone executable, run:

.. code-block:: bash

    mkdir build
    cd build
    cmake -DBUILD_PYTHON=Off -DBUILD_EXECUTABLE=On ..
    make
    make install

This will build the Fortran library and install it in the default location (``~/.local/raffle``).

The default installation is the parallel version of RAFFLE.
To install the serial version of RAFFLE, run the following commands:

.. code-block:: bash

    mkdir build
    cd build
    cmake -DBUILD_PYTHON=Off -DBUILD_EXECUTABLE=On -DCMAKE_BUILD_TYPE=Serial ..
    make
    make install


Installing on MacOS (Homebrew)
==============================

RAFFLE is developed on Linux and MacOS, and should work on both.
However, there are likely some additional steps required to install RAFFLE on MacOS.
This is because **it is not recommended to rely on the Mac system Python, or Fortran and C compilers**.

The recommended way to install Python, gfortran and gcc on MacOS is to use `Homebrew <https://brew.sh>`_.
First, install Homebrew by following the guide on their website.

Once Homebrew is installed, you can install the required dependencies by running:

.. code-block:: bash

    brew install python
    brew install gcc
    brew install gfortran
    export CC=$(brew --prefix gfortran)
    export FC=$(brew --prefix gcc)

Confirm a successful Python installation by running:

.. code-block:: bash

    python --version
    whereis python

This should show the correct Python version (3.11 or later) and path.

Next, if you are using ``pip``, then the following command is found to result in the least issues:

.. code-block:: bash

    python -m pip install --upgrade .[ase]

This ensures that the correct Python version is being called, and that the correct version of ``pip`` is being used.



Testing the installation
=========================

Currently, installation testing is only available for the Fortran library and executable.

Both methods below run the same set of tests (found in the ``tests`` directory of the repository), and should give the same results.
These unit tests are designed to test the core functionality of the Fortran code and ensure that it is working as expected.

.. note::
    The Python package does not currently have a test suite, but this is planned for a future release.
    The functionality of the Python package is provided by the Fortran code.
    The Python package is just a wrapper around the Fortran code, so if the Fortran code is working, then the Python package should also work.
    However, these wrapper functions do definitely need a test suite, and this is planned for a future release.

Testing with fpm
----------------

To test the installation of the Fortran library, navigate to the root directory of the repository and run:

.. code-block:: bash

    fpm test

This will run the test suite for the Fortran library.

Testing with cmake
------------------

To test the installation of the Fortran library, navigate to the directory where the library was built.
If the installation instructions above were followed, this will be in the ``build`` directory within the repository.

To run the test suite, run:

.. code-block:: bash

    ctest

Testing with pytest
-------------------

To test the installation of the Python library, navigate to the repository root directory and run:

.. code-block:: bash

    pytest

This will run the unit tests for the Python wrapper, as well as compile and run the Fortran unit tests.
