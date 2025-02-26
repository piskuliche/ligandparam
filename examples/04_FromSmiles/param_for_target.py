import os
import logging
import sys

from pathlib import Path
import shutil as sh

from ligandparam.io.smiles  import *
from ligandparam.recipes import BuildLigand

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

# Set the leaprc used for the calculation
leaprc = []
leaprc.append("leaprc.RNA.OL3")
leaprc.append("leaprc.gaff2")
leaprc.append("leaprc.water.tip4pew")


reference_structure = "1y27.pdb"
reference_resname = "GUN"

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
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setLevel(logging.INFO)
logger.addHandler(stream_handler)

# baseoptions = {
#     "name": None,
#     "nproc": 12,
#     "mem": "60GB",
#     "net_charge": 0,
#     "atom_type": "gaff2",
#     "leaprc": leaprc,
#     "target_pdb": reference_structure,
#     "force_gaussian_rerun": False
# }

for resname, mol in example_set.items():
    # Generate the PDB from SMILES
    pdb = PDBFromSMILES(resname, mol)
    pdb.mol_from_smiles()
    pdb.write_pdb(f"{resname}_input.pdb")
    
    new = RenamePDBTypes(f"{resname}_input.pdb", mol)
    new.add_mol("1y27_lig.pdb")
    new.rename_by_reference()
     
    # Make a directory for the molecule and cd into it.
    newdir = Path(f"{mol}")
    newdir.mkdir(exist_ok=True)
    
    sh.copyfile(reference_structure, f"{mol}/{reference_structure}")
    sh.copyfile(f"{mol}.pdb", f"{mol}/{mol}.pdb")
    os.system(f"sed -i -e 's@{reference_resname}@{mol}@g' {mol}/{reference_structure}")

    os.chdir(newdir)
    # Do the build
    baseoptions["name"] = mol
    build = BuildLigand(inputoptions=baseoptions)
    build.setup()
    build.list_stages()
    #build.execute(dry_run=False)

    os.chdir(Path(".."))
