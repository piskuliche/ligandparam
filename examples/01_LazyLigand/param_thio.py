#!/usr/bin/env python

# Import the module
from ligandparam.recipes import LazyLigand

inputoptions = {
    'base_name': 'thiophenol',
    'netcharge': 0,
    'mem': '60GB',
    'nproc': 12
}

# Load the pdb as a instance of the LazyLigand class
test = LazyLigand(inputoptions=inputoptions)

# Select the pre-initialized stages for Lazy Ligand
test.setup()

# List the stages out to the user
test.list_stages()

# Execute the stages in order.
test.execute(dry_run=False)
