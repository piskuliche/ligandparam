import shutil
import os

from pathlib import Path

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

for i, molec in enumerate(example_set):
    # Generate the PDB from SMILES
    pdb = PDBFromSMILES(molec, example_set[molec])
    pdb.mol_from_smiles()
    pdb.write_pdb(f"{molec}_input.pdb")
    
    new = RenamePDBTypes(f"{molec}_input.pdb", molec)
    new.add_mol("1y27_lig.pdb")
    new.rename_by_reference()
     
    # Make a directory for the molecule and cd into it.
    newdir = Path(f"{molec}")
    newdir.mkdir(exist_ok=True)
    
    shutil.copyfile(reference_structure, f"{molec}/{reference_structure}")
    shutil.copyfile(f"{molec}.pdb", f"{molec}/{molec}.pdb")
    os.system(f"sed -i -e 's@{reference_resname}@{molec}@g' {molec}/{reference_structure}")

    os.chdir(newdir) 
    # Do the build
    build = BuildLigand(f"{molec}.pdb", netcharge=0, nproc=12,
                         mem='60GB', target_pdb=reference_structure, leaprc=leaprc)
    build.setup()
    build.list_stages()
    #build.execute(dry_run=False)

    os.chdir(Path(".."))
