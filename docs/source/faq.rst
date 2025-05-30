.. _faq:

==========================
Frequently Asked Questions
==========================


Issues with installation
========================

We have found that there are some issues with the installation of RAFFLE on some systems.
These often relate to the Fortran and C compilers.

Pip pointing to incorrect Fortran and C compilers
-------------------------------------------------

If you are using pip to install RAFFLE, you may encounter an issue where pip points to the wrong Fortran and C compilers (this can happen with both local installation and installation from PyPI).
This can happen if you have multiple versions of gfortran and gcc installed on your system, or if you have installed them using a package manager such as Homebrew.
To fix this, you can set the environment variables ``CC`` and ``FC`` to point to the correct compilers before running pip.
To do this, run the following commands:

.. code-block:: bash

    export CC=$(which gfortran)
    export FC=$(which gfortran)
    export CXX=$(which g++)

This will set the environment variables to point to the correct compilers (the ones you currently have in your path).

However, we have found this to not always work either.
So, another option is to run the following line to install RAFFLE:

.. code-block:: bash

    CC=<PATH_TO_C_COMPILER> FC=<PATH_TO_FORTRAN_COMPILER> CXX=<PATH_TO_CXX_COMPILER> pip install --upgrade raffle

Installing on MacOS (Homebrew)
------------------------------

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

    python -m pip install --upgrade .

This ensures that the correct Python version is being called, and that the correct version of ``pip`` is being used.


Issues with running RAFFLE
==========================

Sometimes, RAFFLE installs correctly, but you may encounter issues when running it.
This is often due to missing dependencies or incorrect environment variables.

Missing Fortran libraries
--------------------------

Sometimes, installation works fine, but you may encounter an error when running RAFFLE that looks like this:

.. code-block:: text

    ImportError: libgfortran.so.5: cannot open shared object file: No such file or directory

To fix this, you need to install the correct version of the gfortran library.
If you are using a conda environment, you can do this by running the following command:

.. code-block:: bash

    conda install -c conda-forge libgfortran-ng libgfortran5


General
=======

.. _cite:

How to cite RAFFLE?
--------------------

If you use RAFFLE in your research, please cite the following papers:


    | Joe Pitfield, Ned Thaddeus Taylor, Steven Paul Hepplestone,
    | Predicting Phase Stability at Interfaces,
    | Phys. Rev. Lett. 132, 066201, 2024.
    | doi: 10.1103/PhysRevLett.132.066201

and

    | Ned Thaddeus Taylor, Joe Pitfield, Steven Paul Hepplestone,
    | RAFFLE: Active learning accelerated interface structure prediction,
    | arXiv, 2504.02528v1, 2025.
    | doi: 10.48550/arXiv.2504.02528

BibTex (:git:`bibliography <docs/RAFFLE.bib>`):

.. literalinclude:: ../RAFFLE.bib
