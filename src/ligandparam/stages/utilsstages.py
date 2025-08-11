from rdkit import Chem

__all__ = ("set_atom_pdb_info", )

def set_atom_pdb_info(mol: Chem.Mol, resname: str = "LIG") -> Chem.Mol:
    """
    Sets PDB residue information for all atoms in a molecule.

    Parameters
    ----------
    mol : Chem.Mol
        The RDKit molecule to modify.
    resname : str, optional
        The residue name to set (default is "LIG").

    Returns
    -------
    Chem.Mol
        The modified molecule with PDB residue information set.
    """

    mol.SetProp("_Name", resname)
    mi = Chem.AtomPDBResidueInfo()
    mi.SetResidueName(resname)
    mi.SetResidueNumber(1)
    mi.SetOccupancy(0.0)
    mi.SetTempFactor(0.0)
    [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
    return mol