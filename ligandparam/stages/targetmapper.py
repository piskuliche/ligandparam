import MDAnalysis as mda

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

from ligandparam.stages.abstractstage import AbstractStage


class StageMapTarget(AbstractStage):
    def __init__(self, name, inputoptions=None) -> None:
        """ This class is used to run a basic leap calculations on the ligand.
        
        Parameters
        ----------
        name : str
            The name of the stage
        inputoptions : dict
            The input options
        """

        self.name = name
        self._parse_inputoptions(inputoptions)
        self._add_required(f"{self.base_name}.frcmod")
        self._add_required(f"{self.base_name}.resp.mol2")
        self._add_required(f"{self.target_pdb}")
        self._add_output(f"{self.base_name}.pdb")

        return
    

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage
    

    def _execute(self, dry_run=False):
        """ Map the ligand into the target pdb file using an alignment with the MCS method

        This stage will take the ligand and the target pdb file, and align the ligand to the target pdb file using the MCS method. The ligand-target 
        complex will be saved as a pdb file named after the ligand.
         
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        
        """
        pdb = Chem.MolFromPDBFile(self.target_pdb, removeHs=True)
        lig_mol = Chem.MolFromPDBFile(f"{self.pdb_filename}")

        # Extract mol from pdb (ligand)
        mol = Chem.MolFromPDBBlock("")
        for res in Chem.GetMolFrags(pdb, asMols=True):
            if res.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName() == "LIG":
                mol_from_pdb = res

        # Find the MCS
        mcs = rdFMCS.FindMCS([lig_mol, mol_from_pdb])
        common_smarts = mcs.smartsString
        common_mol = Chem.MolFromSmarts(common_smarts)

        # Get the atom indices of the matches
        match_lig = lig_mol.GetSubstructMatch(common_mol)
        match_pdb = mol_from_pdb.GetSubstructMatch(common_mol)

        # Align the mol
        AllChem.AlignMol(lig_mol, mol_from_pdb, atomMap=list(zip(match_lig, match_pdb)))
        
        Chem.MolToPDBFile(lig_mol, f"{self.base_name}.pdb")

        return
    
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")