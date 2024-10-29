Installation Instructions
=========================

To install the ligand-param package, currently you must use git to clone the repository and install the package. 
This will be updated in the future to allow for installation via pip/anaconda.

Git Instructions
----------------

Packaging for the parametrization script will be improved at a later date, but in the short term,
you can install the package by cloning the repository and installing it using pip into a conda environment with Python 3.10.

.. code-block:: bash

    git clone https://github.com/piskuliche/ligandparam.git
    cd ligandparam
    conda -create -n ligandparam python=3.10
    conda install -c conda-forge parmed
    pip install .

.. note:: Right now the package is only compatible with Python 3.10. Parmed is required for the package to work; however, their versioning is invalid for newer versions of Python.