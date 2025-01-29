#!/usr/bin/env python
from pathlib import Path

# Import the module
from ligandparam.recipes import LazyLigand

inputoptions = {
    'net_charge': 0,
    'mem': '60GB',
    'nproc': 12,
    "atom_type": "gaff2",
    "gaussian_root": "/home/pb777/GAUSSIAN",
    "gauss_exedir": "/home/pb777/GAUSSIAN/g16/bsd:/home/pb777/GAUSSIAN/g16",
    "gaussian_binary": "/home/pb777/GAUSSIAN/g16/g16",
    "gaussian_scratch": "/home/pb777/GAUSSIAN/g16/scratch",

}

cwd = Path("./examples/01_LazyLigand/")
# Load the pdb as a instance of the LazyLigand class
test = LazyLigand(name=cwd / "thiophenol.pdb", cwd=cwd, **inputoptions) -

# Select the pre-initialized stages for Lazy Ligand
test.setup()

# List the stages out to the user
test.list_stages()

# Execute the stages in order.
test.execute(dry_run=False)
