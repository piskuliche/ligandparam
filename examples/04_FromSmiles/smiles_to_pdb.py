import logging
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from ligandparam import __version__
from ligandparam.recipes import LazierLigand
from ligandparam.stages import StageSmilesToPDB

"""
In this example, we're going to parametrize 6 ligands in parallel, starting from their SMILES strings.
"""

def set_file_logger(
        logfilename: Path, logname: str = None, filemode: str = "a"
) -> logging.Logger:
    if logname is None:
        logname = Path(logfilename).stem
    logger = logging.getLogger(logname)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "{asctime} - {levelname} - {version} {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
        defaults={"version": __version__},
    )
    file_handler = logging.FileHandler(filename=logfilename, mode=filemode)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def worker(mol:str, resname:str, data_cwd: Path):
    binder_dir = data_cwd / resname
    binder_dir.mkdir(exist_ok=True)
    binder_pdb = binder_dir / f"{resname}.pdb"
    logger = set_file_logger(binder_dir / f"{resname}.log", filemode="w")
    prestage = StageSmilesToPDB(f"build_{resname}", mol, binder_dir, out_pdb=binder_pdb, resname=resname, logger=logger)

    recipe = LazierLigand(
        in_filename=binder_pdb,
        cwd=binder_dir,
        atom_type="gaff2",
        charge_model = "bcc",
        logger=logger,
        molname=resname,
    )
    recipe.setup()
    recipe.insert_stage(prestage, "Initialize")
    recipe.execute()

    return binder_pdb

# Here is an initial set of molecules
example_set = {
    "F3G": "O=C1NC(C(F)(F)F)=NC2=C1N=CN2",
    "NOG": "O=C1NC(NOC)=NC2=C1N=CN2",
    "DOG": "O=C1NC(NC(C)=O)=NC2=C1N=CN2",
    "NNG": "O=C1NC(NNC)=NC2=C1N=CN2",
    "ORG": "O=C1NC(NOCC2=CC=CC=C2)=NC3=C1N=CN3",
    "LIG": "O=C1NC(NNC2=CC=CC=C2)=NC3=C1N=CN3"
}


cwd = Path(".").resolve()
with ProcessPoolExecutor(max_workers=6) as ex:
    futuros = {}
    for resname, mol in example_set.items():
        futu = ex.submit(worker, mol, resname, cwd)
        # Make sure `resname`s are unique!
        futuros[resname] = futu

    for node, futu in futuros.items():
        try:
            futu.result()
            print(f"Got {node}")
        except Exception as e:
            print(f"Failed {node}: {e}")
            continue
