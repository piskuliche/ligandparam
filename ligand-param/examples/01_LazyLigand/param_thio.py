#!/usr/bin/env python

# Import the module
from ligand_param.recipes import LazyLigand

# Load the pdb as a instance of the LazyLigand class
test = LazyLigand('thiophenol.pdb', netcharge=0,nproc=12,mem='60GB')

# Select the pre-initialized stages for Lazy Ligand
test.setup()

# List the stages out to the user
test.list_stages()

# Execute the stages in order.
test.execute(dry_run=False)
