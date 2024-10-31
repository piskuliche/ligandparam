#!/usr/bin/env python

# Import the module
from ligandparam.recipes import FreeLigand

inputoptions = {
    'base_name': 'thiophenol',
    'netcharge': 0,
    'mem': '60GB',
    'nproc': 12
}

# Load the pdb as a instance of the FreeLigand class
test = FreeLigand('thiophenol.pdb', netcharge=0,nproc=12,mem='60GB')

# Select the pre-initialized stages for Lazy Ligand
test.setup()

# List the stages out to the user

test.list_stages()


test.remove_stage("Normalize")

from ligandparam.stages import StageNormalizeCharge

test.add_stage(StageNormalizeCharge("Normalize2", inputoptions=inputoptions, orig_mol2=test.base_name+".resp.mol2",
						  new_mol2=test.base_name+".resp.mol2"))

test.insert_stage(StageNormalizeCharge("Normalize3", inputoptions=inputoptions, orig_mol2=test.base_name+".resp.mol2",
                                                  new_mol2=test.base_name+".resp.mol2"),"Normalize2")

test.execute()

