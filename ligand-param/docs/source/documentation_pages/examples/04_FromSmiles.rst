Example 04: Building Systems from SMILES
========================================

This example demonstrates how to use the ligand_param package to generate pdb files from SMILES strings,
and use them to parametrize ligands. This example will use the FreeLigand class to parametrize a series of ligands
from a SMILES strings. This example will also demonstrate how to BUILD these ligands into starting
structures for Molecular Dyanmics simulations.

Learning Outcomes:
------------------

1) Learn how to generate pdb files from SMILES strings, and use them to parametrize ligands.
2) Demonstrate automation of batches of ligand parametrizations using python scripting.
3) Demonstrate how to build parm7 and rst7 files for Molecular Dynamics simulations.

Files 
-----
The files for this example can be found in the `LigandParameterization/ligand_param/examples/04_FreeLigand` directory of the source code.


Tutorial 
--------

As with previous examples, we need to import the necessary modules and classes from the ligand_param package. 

.. code-block:: python

    import shutil
    import os

    from pathlib import Path

    from ligand_param.io.smiles  import *
    from ligand_param.recipes import BuildLigand

Here, we are additionally using the pathlib libraries to manage file paths, and the shutil library to copy files.

Next, we define a set of moleucles and their corresponding SMILES strings to a dictionary called example_set. Likewise, 
we define the force-field parameters to be used in the calculation by adding them to a leaprc list. We also define the reference 
structure and residue name IN THE TARGET PDB.

.. code-block:: python
    
    # Set the leaprc used for the calculation
    leaprc = []
    leaprc.append("leaprc.RNA.OL3")
    leaprc.append("leaprc.gaff2")
    leaprc.append("leaprc.water.tip4pew")


    reference_structure = "1y27.pdb"
    reference_resname = "GUN"


Now we loop over the example_set dictionary and do the calculation. 

.. code-block:: python

    for i, molec in enumerate(example_set):
        # (1) Generate the PDB from SMILES
        pdb = PDBFromSMILES(molec, example_set[molec])
        pdb.mol_from_smiles()
        pdb.write_pdb(f"{molec}_input.pdb")
        
        # (2) Rename the atom types to match the reference structure
        new = RenamePDBTypes(f"{molec}_input.pdb", molec)
        new.add_mol("1y27_lig.pdb")
        new.rename_by_reference()
        
        # (3) Make a directory for the molecule and cd into it.
        newdir = Path(f"{molec}")
        newdir.mkdir(exist_ok=True)
        shutil.copyfile(reference_structure, f"{molec}/{reference_structure}")
        shutil.copyfile(f"{molec}.pdb", f"{molec}/{molec}.pdb")
        os.system(f"sed -i -e 's@{reference_resname}@{molec}@g' {molec}/{reference_structure}")

        #(4) Change directory to the new directory and build the systems
        os.chdir(newdir) 
        # Do the build
        build = BuildLigand(f"{molec}.pdb", netcharge=0, nproc=12,
                            mem='60GB', target_pdb=reference_structure, leaprc=leaprc)
        build.setup()
        build.list_stages()

        #build.execute(dry_run=False)

        os.chdir(Path(".."))

Here there are a few things to point out. 

    1) We generate the pdb file from the SMILES string using the PDBFromSMILES class. 
    2) We rename the atom types in the pdb file to match the reference structure using the RenamePDBTypes class. So if one of your atoms in the reference structure is atom N9, then the matching atom in the input pdb file will be renamed to N9. This will ensure that the atom types are shared for common substructures.
    3) We make a directory and copy files into it
    4) We change directory to the new directory and build the ligand using the BuildLigand class.

The BuildLigand class is similar to the FreeLigand class, but it is used to build ligands from pdb files. 
The setup method is used to initialize the stages, and the list_stages method is used to list the stages. 
The execute method is used to run the stages in order.

The main difference is that this class calls a new stage called :class:`ligand_param.recipes.StageBuild`, which builds the ligand into a 
starting structure for MD simulations. 

This class adds the following stages to the pipeline within its setup_method. 

.. code-block:: python

    StageBuild("BuildGas", base_cls=self, build_type='gas'),
    StageBuild("BuildAq", base_cls=self, build_type='aq', concentration=0.14),
    StageBuild("BuildTarget", base_cls=self, build_type='target', target_pdb=self.target_pdb)

The first of these, builds just a gas-phase parm7 and rst7 file. The second builds a parm7 and rst7 file for the ligand in water with a salt concentration
of 0.14M. The third builds a parm7 and rst7 file for the ligand in the target pdb file (aka a protein ligand system or a protein rna system).


Full code
---------

.. code-block:: python

    import shutil
    import os

    from pathlib import Path

    from ligand_param.io.smiles  import *
    from ligand_param.recipes import BuildLigand


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
        # (1) Generate the PDB from SMILES
        pdb = PDBFromSMILES(molec, example_set[molec])
        pdb.mol_from_smiles()
        pdb.write_pdb(f"{molec}_input.pdb")
        
        # (2) Rename the atom types to match the reference structure
        new = RenamePDBTypes(f"{molec}_input.pdb", molec)
        new.add_mol("1y27_lig.pdb")
        new.rename_by_reference()
        
        # (3) Make a directory for the molecule and cd into it.
        newdir = Path(f"{molec}")
        newdir.mkdir(exist_ok=True)
        shutil.copyfile(reference_structure, f"{molec}/{reference_structure}")
        shutil.copyfile(f"{molec}.pdb", f"{molec}/{molec}.pdb")
        os.system(f"sed -i -e 's@{reference_resname}@{molec}@g' {molec}/{reference_structure}")

        #(5) Change directory to the new directory
        os.chdir(newdir) 
        # Do the build
        build = BuildLigand(f"{molec}.pdb", netcharge=0, nproc=12,
                            mem='60GB', target_pdb=reference_structure, leaprc=leaprc)
        build.setup()
        build.list_stages()

        # (6) Execute the build
        #build.execute(dry_run=False)

        os.chdir(Path(".."))
