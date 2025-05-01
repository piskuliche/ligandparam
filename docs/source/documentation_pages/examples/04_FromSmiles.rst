Example 04: Building Systems from SMILES
========================================

This example demonstrates how to use the ligandparam package to generate pdb files from SMILES strings,
and use them to parametrize ligands. This example will use the FreeLigand class to parametrize a series of ligands
from a SMILES strings in parallel. This example will also demonstrate how to BUILD these ligands into starting
structures for Molecular Dyanmics simulations.

Learning Outcomes:
------------------

1) Learn how to generate pdb files from SMILES strings, and use them to parametrize ligands.
2) Demonstrate automation of batches of ligand parametrizations using python scripting.
3) Demonstrate how to build parm7 and rst7 files for Molecular Dynamics simulations.

Files 
-----
The files for this example can be found in the `LigandParameterization/examples/04_FromSmiles` directory of the source code.


Tutorial 
--------

As with previous examples, we need to import the necessary modules and classes from the ligandparam package. 

.. literalinclude :: ../../../../examples/04_FromSmiles/smiles_to_pdb.py
    :language: python
    :end-at: from ligandparam.stages

To work in parallel, we will define a worker and a logger for each task. 

.. literalinclude :: ../../../../examples/04_FromSmiles/smiles_to_pdb.py
    :language: python
    :start-at: set_file_logger
    :end-at: Here is an initial

Next, we define a set of moleucles and their corresponding SMILES strings to a dictionary called example_set. Likewise, 
we define the force-field parameters to be used in the calculation by adding them to a leaprc list. We also define the reference 
structure and residue name IN THE TARGET PDB.


.. literalinclude :: ../../../../examples/04_FromSmiles/smiles_to_pdb.py
    :language: python
    :start-at: example_set = {
    :end-at: "LIG":


Lastly, we run the calculations.

.. literalinclude :: ../../../../examples/04_FromSmiles/smiles_to_pdb.py
    :language: python
    :start-at: cwd = Path


Full code
---------


.. literalinclude :: ../../../../examples/04_FromSmiles/smiles_to_pdb.py
    :language: python

