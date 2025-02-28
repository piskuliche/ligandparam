import logging
import sys

from pathlib import Path

from ligandparam import LazierLigand
from ligandparam.stages.stagesmiles import StageSmilesToPDB

# Here is an initial set of molecules 
example_set = {
    "F3G": "O=C1NC(C(F)(F)F)=NC2=C1N=CN2",
    "NOG": "O=C1NC(NOC)=NC2=C1N=CN2",
    "DOG": "O=C1NC(NC(C)=O)=NC2=C1N=CN2",
    "NNG": "O=C1NC(NNC)=NC2=C1N=CN2",
    "ORG": "O=C1NC(NOCC2=CC=CC=C2)=NC3=C1N=CN3",
    "NNG": "O=C1NC(NNC)=NC2=C1N=CN2",
    "LIG": "O=C1NC(NNC2=CC=CC=C2)=NC3=C1N=CN3"
}

cwd = Path(".").resolve()
# Send output to stdout
logger = logging.getLogger("mylog")
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setLevel(logging.INFO)
logger.addHandler(stream_handler)

recipes = []
# Generate the PDB from the SMILES
for resname, mol in example_set.items():
    pdb = cwd / f"{resname}.pdb"
    prestage = StageSmilesToPDB(f"build_{resname}", mol, cwd, out_pdb=pdb, resname=resname, logger=logger)
    recipe = LazierLigand(
        in_filename=pdb,
        cwd=cwd,
        atom_type="gaff2",
        logger=logger,
        molname=resname,
    )
    recipe.setup()
    recipe.insert_stage(prestage, "Initialize")
    recipes.append(recipe)

for r in recipes:
    r.execute()