Installation Instructions
=========================

To install the ligand-param package, currently you must use git to clone the repository and install the package. 
This will be updated in the future to allow for installation via pip/anaconda.

Git Instructions
----------------

To install the ligandparam package from GitHub, we use miniforge (a minimal version of conda) to create a new environment and install the package using
the pre-built conda pacakge described in `env.yaml`. 

For this, to install miniforge you need to follow the instructions in the `miniforge` github page, but minimally this involves
downloading the installer script and running it. Miniforge is located here: `Miniforge <https://github.com/conda-forge/miniforge?tab=readme-ov-file>`_.

A simple install for miniforge is as follows:

.. code-block:: bash
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash Miniforge3-Linux-x86_64.sh

Once miniforge is installed, you can create a new environment and install the ligandparam package using the following commands:

.. code-block:: bash

    source ~/.bashrc # to update environment if you just installed miniforge.
    git clone https://github.com/piskuliche/ligandparam.git
    cd ligandparam
    mamba env create -f env.yaml
    conda activate ligandparam
    pip install -e .


This will install the ligandparam package in the current environment, making it available for use in python scripts.

As of now, this package also requires antechamber, tleap, and g16 to be accessible within your PATH environment. 


    
