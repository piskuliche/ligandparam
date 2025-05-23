import sys
from pathlib import Path

# Import the module
from ligandparam.recipes import LazyLigand
from ligandparam.stages import StageNormalizeCharge

# Environment variables for Gaussian. If your environment is already set up, you can ignore this.
gaussian_paths = {
    "gaussian_root": "/home/.../GAUSSIAN",
    "gauss_exedir": "/home/.../GAUSSIAN/g16/bsd:/home/.../GAUSSIAN/g16",
    "gaussian_binary": "/home/.../GAUSSIAN/g16/g16",
    "gaussian_scratch": "/home/.../GAUSSIAN/g16/scratch",
}

cwd = Path(sys.argv[0]).resolve().parent

# Load the pdb as a instance of the FreeLigand class
test = LazyLigand(
    in_filename=cwd / "thiophenol.pdb",
    cwd=cwd,
    net_charge=0,
    atom_type="gaff2",
    logger="stream",
    # antechamber will name your residue 'MOL' by default, and we follow that standard by default,
    # so you probably want to set it yourself:
    molname="LIG",
    # **gaussian_paths,
)

# Select the pre-initialized stages for Lazy Ligand
test.setup()

# List the stages out to the user
test.list_stages()

test.remove_stage("Normalize1")

test.insert_stage(
    StageNormalizeCharge("mynormalization", main_input="thiophenol.initial.mol2", cwd=cwd, net_charge=0,
                         out_mol2="thiophenol.initial.mol2"),
    "MinimizeLowTheory",
)

# List the stages after inserting.

test.list_stages()

test.execute(nproc=12, mem="8")
