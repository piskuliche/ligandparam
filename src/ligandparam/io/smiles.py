import rdkit
import re
import shutil
import MDAnalysis as mda
import numpy as np

from rdkit.Chem import rdFMCS

class PDBFromSMILES:
    """ This class generates a PDB file from a SMILES string. 
    
    Parameters:
    ----------
    resname : str
        The residue name of the molecule.
    smiles : str
        The SMILES string of the molecule.
    

    Attributes:
    ----------
    resname : str
        The residue name of the molecule.
    smiles : str
        The SMILES string of the molecule.
    mol : rdkit.Chem.Mol
        The RDKit molecule object.
    pdb_filename : str
        The name of the PDB file.
    """
    def __init__(self, resname, smiles):
        self.resname = resname
        self.smiles = smiles
        self.mol = None
        return

    
    def write_pdb(self, filename, randomSeed=0xf00d):
        """ Write the PDB file using rdkit to filename. 
        
        Parameters:
        ----------
        filename : str
            The name of the PDB file.
        randomSeed : hex
            The random seed for the embedding algorithm.
        
        Returns:
        --------
        None
        """
        params = rdkit.Chem.AllChem.ETKDGv3()
        params.randomSeed = randomSeed
        rdkit.Chem.AllChem.EmbedMolecule(self.mol, params)
        rdkit.Chem.rdmolfiles.MolToPDBFile(self.mol, filename)
        self.pdb_filename = filename
        clean_pdb(self.pdb_filename, self.resname)
        return
    
    def mol_from_smiles(self, addHs=True):
        """ Generate a molecule from a SMILES string. 
        
        Parameters:
        -----------
        addHs : bool
            Whether to add hydrogens to the molecule.
            
        Returns:
        -------
        None
        """
        mol = rdkit.Chem.MolFromSmiles(self.smiles)
        if addHs:
            mol = rdkit.Chem.rdmolops.AddHs(mol)
        self.mol = mol
        return 
    
class MolFromPDB:
    def __init__(self, pdb_filename, removeHs=False):
        self.remove_Hs = removeHs
        self.pdb_filename = pdb_filename
        self._rdkit_representation()
        self._mda_representation()
        return

    
    def _rdkit_representation(self):
        """ Generate an RDKit molecule from a PDB file. """
        self.rdkit_mol = rdkit.Chem.rdmolfiles.MolFromPDBFile(self.pdb_filename, removeHs=self.remove_Hs)
        return
    
    def _mda_representation(self):
        """ Generate an MDAnalysis Universe from a PDB file. """
        self.mda_universe = mda.Universe(self.pdb_filename)
        return
    
    def resname(self):
        """ Get the residue name from the PDB file. """
        return self.mda_universe.atoms.residues.resnames[0]
    
    def names(self):
        return self.mda_universe.atoms.names
    
    def elements(self):
        return self.mda_universe.atoms.elements
    
    def write_pdb(self, filename):
        """ Write the PDB file using MDAnalysis to filename.
        
        Parameters:
        ----------
        filename : str
            The name of the PDB file.
        
        """
        self.mda_universe.atoms.write(filename)
        clean_pdb(filename, self.resname())
        return
    
    
class RenamePDBTypes:
    """ This class renames the atom types in a PDB file based on a reference PDB file.
    
    Parameters:
    ----------
    primary_pdb : str
        The name of the primary PDB file.
    resname : str
        The residue name of the molecule.
    
    Attributes:
    -----------
    primary_pdb : str
        The name of the primary PDB file.
    resname : str
        The residue name of the molecule.
    mols : list
        A list of MolFromPDB objects.
    mcs_mol : rdkit.Chem.Mol
        The RDKit molecule object of the common substructure.
    
    
    """
    def __init__(self, primary_pdb, resname):
        self.primary_pdb = primary_pdb
        self.mols = []
        self.mols.append(MolFromPDB(primary_pdb, removeHs=False))
        self.resname = resname
        return
    
    def add_mol(self, mol_pdb):
        """ Add a molecule to the RenamePDBTypes object. """
        self.mols.append(MolFromPDB(mol_pdb))
        return
    
    def rename_by_reference(self):
        """ Rename the PDB files by the reference PDB file. """
        if len(self.mols) != 2:
            raise ValueError("ERROR: Only two molecules can be compared for reference.")
        st_comm, rf_comm = self.common_atoms()
        codes = [f"C{i}" for i in range(1, 41)]
        available_names = set(self.mols[0].names()) | set(self.mols[1].names()) | set(codes)
        new_names = np.zeros_like(self.mols[0].names())
        print(self.mols[0].names()[st_comm])
        print(self.mols[1].names()[rf_comm])
        for i, j in zip(st_comm, rf_comm):
            new_names[i] = self.mols[1].names()[j]
            available_names.remove(self.mols[1].names()[j])
        for atom in self.mols[0].mda_universe.atoms:
            if atom.index not in st_comm:
                renamed=False
                for key in available_names:
                    cleaned_key = self._split_letters_numbers(key)
                    if self.mols[0].elements()[atom.index] == cleaned_key[0]:
                        new_names[atom.index] = key
                        available_names.remove(key)
                        renamed=True
                        break
                if not renamed:
                    raise ValueError("ERROR: Could not rename atom. ")
        self.mols[0].mda_universe.atoms.names = new_names
        shutil.copyfile(self.mols[0].pdb_filename, 'original_'+self.mols[0].pdb_filename)
        self.mols[0].write_pdb(f"{self.mols[0].resname()}.pdb")
                    
        return
    
    def find_mcs(self):
        """ Find the common substructure between the molecules. """
        mcs = rdFMCS.FindMCS([mol.rdkit_mol for mol in self.mols])
        self.mcs_mol = rdkit.Chem.rdmolfiles.MolFromSmarts(mcs.smartsString)
        return
    
    def common_atoms(self):
        """ Get the common atoms between the molecules. """
        self.find_mcs()
        return [np.array(mol.rdkit_mol.GetSubstructMatch(self.mcs_mol)) for mol in self.mols]

    def _split_letters_numbers(self, s):
        """Split a string into letters and numbers."""
        match = re.match(r"([a-zA-Z]+)([0-9]+)", s)
        if match:
            return match.groups()
        else:
            return None


def clean_pdb(pdb_filename, resname):
    if len(resname) != 3:
        raise ValueError("Resname must be 3 characters")
    lines = []
    with open(pdb_filename, "r") as f:
        lines = f.readlines()
    with open(pdb_filename, "w") as f:
        for line in lines:
            line = line.replace("UNL", resname)
            line = line.replace("SYST", "    ")
            if line.startswith("ATOM") or line.startswith("HETATM"):
                f.write(line)

    return
    
if __name__ == "__main__":
    # Create the PDBFromSMILES object
    pdb = PDBFromSMILES("F3G", "O=C1NC(C(F)(F)F)=NC2=C1N=CN2")
    
    # Generate the molecule
    mol = pdb.mol_from_smiles()
    
    # Write the PDB file
    pdb.write_pdb(mol, "test.pdb")

    new = RenamePDBTypes("test.pdb", "F3G")
    new.add_mol("align.pdb")
    new.rename_by_reference()
