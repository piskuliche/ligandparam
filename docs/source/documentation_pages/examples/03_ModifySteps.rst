Example 03: Utilizing Step Modifications
========================================

This example demonstrates how to use the FreeLigand class to perform a multi-state RESP calculation, with modifications made 
to the steps in the pipeline. This tutorial shows briefly how to add, remove, and insert stages into the pipeline, starting from
a pre-initialized pipeline.

Learning Outcomes:
------------------

1) Learn how top modify existing pipeline recipes by adding, removing, and inserting stages.

Files 
-----
The files for this example can be found in the `LigandParameterization/examples/03_ModifySteps` directory of the source code.


Tutorial 
--------

In this example, we will demonstrate how to take an existing recipe, such as ModifySteps, and then modify the stages in the pipeline. This can 
also be used to create a fully new pipeline by subclassing the Parameterization class and adding your own stages.

As always, the first step is to import the module, load the pdb file into the recipe of your choosing, and use the setup method and list_stages method
to add the stages to the pipeline, and print them to the screen.

.. code-block:: python

    # Import the module
    from ligandparam.recipes import *

    inputoptions = {
        'base_name': 'thiophenol',
        'net_charge': 0,
        'mem': '60GB',
        'nproc': 12
    }
    # Load the pdb as a instance of the FreeLigand class
    test = FreeLigand(inputoptions=inputoptions)

    # Select the pre-initialized stages for Lazy Ligand
    test.setup()

    # List the stages out to the user
    test.list_stages()

Now, unlike before, we are going to remove the Normalize1 stage from the pipeline. If you recall, the Normalize1 stage occurs
after the Initialize stage and before the Minimize stage.

.. code-block:: python

    test.remove_stage("Normalize1")

Running this command removes Normalize1 from the pipeline, and the pipeline will now skip this stage when it is executed (if you use execute).

Stages can also be added to the pipeline, this requires two steps:

1) Import the stage from the module
2) Add the stage to the pipeline

.. code-block:: python

    from ligandparam.stages import StageNormalizeCharge

    test.add_stage(StageNormalizeCharge("Normalize2", inputoptions=inputoptions, 
                                    orig_mol2=test.base_name+".resp.mol2",
						            new_mol2=test.base_name+".resp.mol2"))

Adding stages to the pipeline adds them to the end, so before, the Leap stage was the last stage in the pipeline, but
now the NormalizeA stage is the last stage in the pipeline.

Stages can also be inserted into the pipeline, this requires two steps:

1) Import the stage from the module (if not already imported)
2) Insert the stage into the pipeline

.. code-block:: python

    test.insert_stage(StageNormalizeCharge("Normalize3", inputoptions=inputoptions,
                                        orig_mol2=test.base_name+".resp.mol2",
                                        new_mol2=test.base_name+".resp.mol2"),
                                        "Normalize2")

Here, the Normalize3 stage is inserted into the pipeline before the Normalize2 stage. This is useful if you want to add a stage
to the pipeline, but not at the end of the pipeline.

Finally, the pipeline can be executed using the execute method as before.

.. code-block:: python

    test.execute()

Full code
---------

.. code-block:: python
    
   # Import the module
    from ligandparam.recipes import FreeLigand

    inputoptions = {
        'base_name': 'thiophenol',
        'net_charge': 0,
        'mem': '60GB',
        'nproc': 12
    }

    # Load the pdb as a instance of the FreeLigand class
    test = FreeLigand(inputoptions=inputoptions)

    # Select the pre-initialized stages for Lazy Ligand
    test.setup()

    # List the stages out to the user

    test.list_stages()


    test.remove_stage("Normalize")

    from ligandparam.stages import StageNormalizeCharge

    test.add_stage(StageNormalizeCharge("Normalize2", inputoptions=inputoptions, 
                                        orig_mol2=test.base_name+".resp.mol2",
                                        new_mol2=test.base_name+".resp.mol2"))

    test.insert_stage(StageNormalizeCharge("Normalize3", inputoptions=inputoptions,
                                            orig_mol2=test.base_name+".resp.mol2",
                                            new_mol2=test.base_name+".resp.mol2"),
                                            "Normalize2")

    test.execute()

