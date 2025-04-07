Creating New Recipes
=====================

Recipes are classes that combine the building blocks of the workflow for parametrizing a ligand. 

`LazyLigand`, `FreeLigand`, and `BuildLigand` are examples of recipes that are included in the ligandparam package by default.

However, you might find yourself wanting to generate new recipes that are tailored to your specific needs and/or workflow.

In this example, we will demonstrate how to create a new recipe by subclassing the `Recipe` class and adding your own stages to the pipeline.

Building from Scratch for One Time Use
--------------------------------------

If you are building a recipe for a one-time use, you can use the `Recipe` class directly and add the stages to the pipeline.

Here is an example of a simple recipe that only has two stages, `StageInitialize` and `StageNormalizeCharge`:

.. code-block:: python

    from ligandparam.recipes import Recipe
    from ligandparam.stages import *

    inputoptions = {
        'base_name': 'my_ligand',
        'net_charge': 0,
        'mem': '60GB',
        'nproc': 12
    }

    new_recipe = Recipe(inputoptions=inputoptions)
    new_recipe.add_stage(StageInitialize("Initialize", inputoptions=inputoptions))
    new_recipe.add_stage(StageNormalizeCharge("NormalizeCharge", 
                        inputoptions=inputoptions, 
                        orig_mol2=new_recipe.base_name+".antechamber.mol2",
                        new_mol2=new_recipe.base_name+".antechamber.mol2"))
    new_recipe.execute(dry_run=False)

This script might be useful for generating a mol2 file from a pdb file, creating bcc charges, and then normalizing the charges to zero.


Building from Scratch for Reuse
-------------------------------

If you are building a recipe that you plan to reuse, you should subclass the `Recipe` class and add the stages to the pipeline.

Here is an example of a simple recipe that only has two stages like before, `StageInitialize` and `StageNormalizeCharge`:

.. code-block:: python

    from ligandparam.recipes import Recipe
    from ligandparam.stages import *

    class MyNewRecipe(Recipe):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            return
        def setup(self):
            self.stages = [
                StageInitialize("Initialize", inputoptions=self.inputoptions),
                StageNormalizeCharge("Normalize1", self.inputoptions, 
                                    orig_mol2=self.base_name+".antechamber.mol2", 
                                    new_mol2=self.base_name+".antechamber.mol2")
            ]

Then you can use this recipe just like before.:

.. code-block:: python

    inputoptions = {
        'base_name': 'my_ligand',
        'net_charge': 0,
        'mem': '60GB',
        'nproc': 12
    }

    new_recipe = MyNewRecipe('my_ligand.pdb', netcharge=0, nproc=12, mem='60GB')
    new_recipe.setup()
    new_recipe.execute(dry_run=False)