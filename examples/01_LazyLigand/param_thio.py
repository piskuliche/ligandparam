from pathlib import Path
import sys

from ligandparam.recipes import LazyLigand

cwd = Path(sys.argv[0]).resolve().parent

# Environment variables for Gaussian. If your environment is already set up, you can ignore this.
gaussian_paths = {
    "gaussian_root": "/home/.../GAUSSIAN",
    "gauss_exedir": "/home/.../GAUSSIAN/g16/bsd:/home/.../GAUSSIAN/g16",
    "gaussian_binary": "/home/.../GAUSSIAN/g16/g16",
    "gaussian_scratch": "/home/.../GAUSSIAN/g16/scratch",
}

parametrize_ligand = LazyLigand(
    in_filename=cwd / "thiophenol.pdb",
    cwd=cwd,
    logger="file",
    net_charge=0,
    atom_type="gaff2",
    # antechamber will name your residue 'MOL' by default, and we follow that standard by default,
    # so you probably want to set it yourself:
    molname="LIG",
    # **gaussian_paths,
)

# Set the pre-initialized stages for Lazy Ligand
parametrize_ligand.setup()

# List the stages out to the user
parametrize_ligand.list_stages()

# Run the parametrization
parametrize_ligand.execute(dry_run=False, nproc=12, mem=8)
