Installation Instructions
=========================

To install the ligand-param package, currently you must use git to clone the repository and install the package. 
This will be updated in the future to allow for installation via pip/anaconda.

Git Instructions
----------------

To install via git, clone the repository and install the package.

.. code-block:: bash

    git clone path/to/ligandparam
    cd LigandParamtrize/ligandparam
    pip install .

Note - it is recommended to be using a conda environment with Python 3.10 (for now) to install and run the package.

To create a conda environment, you can use the following command:

.. code-block:: bash

    conda create -n ligandparam python=3.10
    conda activate ligandparam

Then you can install the package as described above. Currently, due to dependencies of the Parmed package, the present package is incompatible
with python versions that use numpy 2.0 or above. 