#!/usr/bin/env python
from pathlib import Path
import logging

# Import the module
from ligandparam.recipes import LazyLigand


# Helper function to set up a file logger
def set_file_logger(logfilename: Path, logname: str = None, filemode: str = 'a') -> logging.Logger:
    if logname is None:
        logname = Path(logfilename).stem
    logger = logging.getLogger(logname)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler = logging.FileHandler(filename=logfilename, mode=filemode)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


cwd = Path(".").resolve()

# Environment variables for Gaussian. If your environment is already set up, you can ignore this.
gaussian_paths = {
    "gaussian_root": "/home/pb777/GAUSSIAN",
    "gauss_exedir": "/home/pb777/GAUSSIAN/g16/bsd:/home/pb777/GAUSSIAN/g16",
    "gaussian_binary": "/home/pb777/GAUSSIAN/g16/g16",
    "gaussian_scratch": "/home/pb777/GAUSSIAN/g16/scratch",
}

logger = set_file_logger(cwd / "lazyligand.log", filemode="w")
parametrize_ligand = LazyLigand(
    in_filename=cwd / "thiophenol.pdb",
    cwd=cwd,
    logger=logger,
    net_charge=0,
    atom_type="gaff2",
    # antechamber will name your residue 'MOL' by default, and we follow that standard by default,
    # so you probably want to set it yourself:
    molname="LIG",
    **gaussian_paths,
)

# Set the pre-initialized stages for Lazy Ligand
parametrize_ligand.setup()

# List the stages out to the user
parametrize_ligand.list_stages()

# Run the parametrization
parametrize_ligand.execute(dry_run=False, nproc=12, mem="8192")
