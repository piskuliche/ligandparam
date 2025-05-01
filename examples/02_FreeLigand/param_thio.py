import sys
from pathlib import Path

from ligandparam.recipes import FreeLigand

cwd = Path(sys.argv[0]).resolve().parent
# Environment variables for Gaussian. If your environment is already set up, you can ignore this.
gaussian_paths = {
    "gaussian_root": "/home/.../GAUSSIAN",
    "gauss_exedir": "/home/.../GAUSSIAN/g16/bsd:/home/.../GAUSSIAN/g16",
    "gaussian_binary": "/home/.../GAUSSIAN/g16/g16",
    "gaussian_scratch": "/home/.../GAUSSIAN/g16/scratch",
}

# Load the pdb as a instance of the FreeLigand class
parametrize_ligand = FreeLigand(
    in_filename=cwd / "thiophenol.pdb",
    cwd=cwd,
    logger="file",
    net_charge=-1,
    atom_type="gaff2",
    molname="LIG",
    # **gaussian_paths,
)

# Select the pre-initialized stages for Lazy Ligand
parametrize_ligand.setup()

# List the stages out to the user
parametrize_ligand.list_stages()

# Execute the stages in order.
parametrize_ligand.execute(dry_run=False, nproc=24, mem=48)
