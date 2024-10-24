.. ligand_param documentation master file, created by
   sphinx-quickstart on Thu Oct 10 16:49:20 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ligand_param's documentation!
========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   module
   stages
   multiresp

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
   from ligand_param.module import *

   # Load the pdb as a instance of the FreeLigand class
   test = FreeLigand('thiophenol.pdb', netcharge=0,nproc=12,mem='60GB')

   # Select the pre-initialized stages for Lazy Ligand
   test.setup()

   # List the stages out to the user
   test.list_stages()

   # Execute the stages in order.
   test.execute(dry_run=False)

The above script will load the pdb file `thiophenol.pdb` and initialize the stages for the ligand parameterization. 

There are a number of different parametrization classes, including: 

- LazyLigand  Which just does a single RESP calculation
- FreeLigand: Which does a multiRESP calculation

These classes are designed to be easily extensible, so you can add your own stages to the pipeline.

These stages are all subclasses of the Parameterization class, which is a subclass of the Driver class.

A key aspect of the code is that you can easily add your own stages to the pipeline, as mentioned above. 

.. code-block:: python
   
   # Import the module
   from ligand_param.module import *

   # Load the pdb as a instance of the FreeLigand class
   test = FreeLigand('thiophenol.pdb', netcharge=0,nproc=12,mem='60GB')

   # Select the pre-initialized stages for Lazy Ligand
   test.setup()

   # List the stages out to the user
   test.list_stages()

   # Remove the Normalize stage
   test.remove_stage("Normalize")

   from ligand_param.stages.charge import StageNormalizeCharge

   # Add a new Normalize stage
   test.add_stage(StageNormalizeCharge("Normalize2", orig_mol2=test.base_name+".resp.mol2",
                     new_mol2=test.base_name+".resp.mol2"))

   # Insert a new Normalize stage before Normalize2
   test.insert_stage(StageNormalizeCharge("Normalize3", orig_mol2=test.base_name+".resp.mol2",
                                                   new_mol2=test.base_name+".resp.mol2"),"Normalize2")

   test.execute()

