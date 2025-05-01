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

.. literalinclude :: ../../../../examples/03_ModifySteps/param_thio.py
    :language: python
    :end-before: test.remove_stage

Now, unlike before, we are going to remove the Normalize1 stage from the pipeline. If you recall, the Normalize1 stage occurs
after the Initialize stage and before the Minimize stage.

.. literalinclude :: ../../../../examples/03_ModifySteps/param_thio.py
    :language: python
    :start-at: test.remove_stage
    :end-at: test.remove_stage


Running this command removes Normalize1 from the pipeline, and the pipeline will now skip this stage when it is executed (if you use execute).

Stages can also be added to the pipeline:

.. literalinclude :: ../../../../examples/03_ModifySteps/param_thio.py
    :language: python
    :start-at: test.insert_stage
    :end-at: "MinimizeLowTheory"

Adding stages to the pipeline adds them before the specified stage, so we added mynormalization before the MinimizeLowTheory stage.


Finally, the pipeline can be executed using the execute method as before.

.. literalinclude :: ../../../../examples/03_ModifySteps/param_thio.py
    :language: python
    :start-after: after inserting

Full code
---------

.. literalinclude :: ../../../../examples/03_ModifySteps/param_thio.py
    :language: python


