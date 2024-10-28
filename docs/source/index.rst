.. ligand_param documentation master file, created by
   sphinx-quickstart on Thu Oct 10 16:49:20 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ligand_param's documentation!
========================================



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Example Script
==============

For detailed examples, please refer to the `examples` directory in the source code.

Here is a simple example of how to use the ligand_param module:

.. code-block:: python 

   #!/usr/bin/env python

   # Import the module
   from ligandparam.recipes import FreeLigand

   # Load the pdb as a instance of the FreeLigand class
   parametrize_ligand = FreeLigand('thiophenol.pdb', netcharge=0,nproc=12,mem='60GB')

   # Select the pre-initialized stages for Lazy Ligand
   parametrize_ligand.setup()

   # List the stages out to the user
   parametrize_ligand.list_stages()

   # Execute the stages in order.
   parametrize_ligand.execute(dry_run=False)

The above script will load the pdb file `thiophenol.pdb` and initialize the stages for the ligand parameterization. 

There are a number of different parametrization classes, including: 

- LazyLigand  Which just does a single RESP calculation
- FreeLigand: Which does a multiRESP calculation
- BuildLigand: Which builds a ligand from a pdb file (and outputs a parm7/rst7 file in gas, aqueous, or a protein/rna target)

These classes are designed to be easily extensible, so you can add your own stages to the pipeline.

These stages are all subclasses of the Parameterization class, which is a subclass of the Driver class.

A key aspect of the code is that you can easily add your own stages to the pipeline, as mentioned above. 

.. code-block:: python
   
   # Import the module
   from ligandparam.recipes import FreeLigand

   # Load the pdb as a instance of the FreeLigand class
   parametrize_ligand = FreeLigand('thiophenol.pdb', netcharge=0,nproc=12,mem='60GB')

   # Select the pre-initialized stages for Lazy Ligand
   parametrize_ligand.setup()

   # List the stages out to the user
   parameterize_ligand.list_stages()

   # Remove the Normalize stage
   parametrize_ligand.remove_stage("Normalize1")

   from ligand_param.stages import StageNormalizeCharge

   # Add a new Normalize stage
   parametrize_ligand.add_stage(StageNormalizeCharge("Normalize2", orig_mol2=test.base_name+".resp.mol2",
                     new_mol2=test.base_name+".resp.mol2"))

   # Insert a new Normalize stage before Normalize2
   parametrize_ligand.insert_stage(StageNormalizeCharge("Normalize3", orig_mol2=test.base_name+".resp.mol2",
                                                   new_mol2=test.base_name+".resp.mol2"),"Normalize2")

   parametrize_ligand.execute()

.. Contents
   ========

.. toctree::
   :maxdepth: 4
   :caption: Documentation
   :numbered:
   :hidden:
   
   ./documentation_pages/overview.rst
   ./documentation_pages/installation.rst
   ./documentation_pages/recipes.rst
   ./documentation_pages/stages.rst
   ./documentation_pages/io.rst
   ./documentation_pages/multiresp.rst
   ./documentation_pages/examples
   