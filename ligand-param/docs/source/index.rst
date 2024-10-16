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

.. code-block:: python 

   from ligand_param.module import *
   from ligand_param.stages.checkcharge import StageCheckCharge

   test = LazyLigand('FBKRP.pdb', netcharge=0,nproc=12,mem='40GB')
   test.setup()
   test.list_stages()
   test.insert_stage(StageCheckCharge("Check1", filename='FBKRP.resp.mol2', filetype='MOL2'),"Leap")
   test.add_stage(StageCheckCharge("Check2", filename='FBKRP.resp.mol2', filetype='MOL2'))
   test.remove_stage("Check1")
   test.execute(dry_run=False)
   test.clean()