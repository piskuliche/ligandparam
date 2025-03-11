import logging
import sys
from pathlib import Path

# Import the module
from ligandparam.recipes import LazyLigand
from ligandparam.stages import StageNormalizeCharge

# Environment variables for Gaussian. If your environment is already set up, you can ignore this.
gaussian_paths = {
    "gaussian_root": "/home/pb777/GAUSSIAN",
    "gauss_exedir": "/home/pb777/GAUSSIAN/g16/bsd:/home/pb777/GAUSSIAN/g16",
    "gaussian_binary": "/home/pb777/GAUSSIAN/g16/g16",
    "gaussian_scratch": "/home/pb777/GAUSSIAN/g16/scratch",
}

cwd = Path(".").resolve()

# Send output to stdout
logger = logging.getLogger("mylog")
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setLevel(logging.INFO)
logger.addHandler(stream_handler)

# Load the pdb as a instance of the FreeLigand class
test = LazyLigand(
    in_filename=cwd / "thiophenol.pdb",
    cwd=cwd,
    net_charge=0,
    atom_type="gaff2",
    logger=logger,
    # antechamber will name your residue 'MOL' by default, and we follow that standard by default,
    # so you probably want to set it yourself:
    molname="LIG",
    **gaussian_paths,
)

# Select the pre-initialized stages for Lazy Ligand
test.setup()

# List the stages out to the user
test.list_stages()

test.remove_stage("Normalize1")

test.insert_stage(
    StageNormalizeCharge(
        "mynormalization",
        cwd=cwd,
        net_charge=0,
        input="thiophenol.initial.mol2",
        out_mol2="thiophenol.initial.mol2",
    ),
    "LazyResp",
)

test.list_stages()

test.execute(nproc=12, mem="8192")
