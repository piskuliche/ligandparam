#!/usr/bin/env python

# Import the module
from ligand_param.module import *

# Load the pdb as a instance of the FreeLigand class
test = FreeLigand('thiophenol.pdb', netcharge=0,nproc=12,mem='60GB')

# Select the pre-initialized stages for Lazy Ligand
test.setup()

# List the stages out to the user
print("The stage system can list current stages with the list_stages() method")
test.list_stages()

print("Steps can be removed using the remove_stage method")
test.remove_stage("Normalize")
print("Note that currently, there is no way of ensuring that")
print("stages that are removed will not break further points")
print("In the script.")

print()
print("Stages can be added just as easily.")
print("Start by importing the stage.")
from ligand_param.stages.charge import StageNormalizeCharge
print("Then add the stage, with any options included using the add_stage option.")
test.add_stage(StageNormalizeCharge("Normalize2", orig_mol2=test.base_name+".resp.mol2",
						  new_mol2=test.base_name+".resp.mol2"))
print("This appends the stage at the END of the list of stages.")
print("Note where this ends up in the list.")

print()
print("You can also insert stages before other stages.")
print("Lets put a stage Normalize3 in front of Normalize2")
print("To do this, lets use the insert_stage method.")
test.insert_stage(StageNormalizeCharge("Normalize3", orig_mol2=test.base_name+".resp.mol2",
                                                  new_mol2=test.base_name+".resp.mol2"),"Normalize2")

print("Note that we are identifying the stage that we want to insert in front of.")

print("Finally, you can call the overall parametrization using")
print("test.execute()")
print("If you want to do this, you can uncomment this in the python script.")
#test.execute()

print("This same procedure can be used with the base Parametrization class")
print("To create a completely customized ligand parametrization, with options")
print("Note - the base Parametrization class does not have a setup method")
print(" so by default, it won't have any stages at the start.")
