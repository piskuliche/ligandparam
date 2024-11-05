import shutil
import os

from pathlib import Path

from ligandparam.io.smiles  import *
from ligandparam.recipes import BuildLigand

# Example default stage list, which could be passed to the disable_stages method to mass remove stages from 
# the recipe. To do that, you would uncomment the line within the for loop marked by a commment.
# Any stage that is set to false will be removed. Note - this is not guaranteed to work, as the stages are
# sometimes dependent on each other.
default_stage_list = {
    "Initialize": True,
    "Normalize1": True,
    "Minimize": True,
    "Rotate": True,
    "GrabGaussianCharge": True,
    "MultiRespFit": True,
    "UpdateCharge": True,
    "Normalize2": True,
    "UpdateNames": True,
    "UpdateTypes": True,
    "ParmChk": True,
    "Leap": True,
    "BuildGas": True,
    "BuildAq": True,
    "BuildTarget": True
}

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


baseoptions = {
    "base_name": None,
    "nproc": 12,
    "mem": "60GB",
    "net_charge": 0,
    "atom_type": "gaff2",
    "leaprc": leaprc,
    "target_pdb": reference_structure,
    "force_gaussian_rerun": False
}

for i, molec in enumerate(example_set):
    # Generate the PDB from SMILES
    pdb = PDBFromSMILES(molec, example_set[molec])
    pdb.mol_from_smiles()
    pdb.draw_mol(f"{molec}_smiles.png")
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
    baseoptions["base_name"] = molec
    build = BuildLigand(inputoptions=baseoptions)

    build.setup()
    # Uncomment the line below to disable all stages
    #build.disable_stages(default_stage_list)

    build.list_stages()
    #build.execute(dry_run=False)

    os.chdir(Path(".."))
