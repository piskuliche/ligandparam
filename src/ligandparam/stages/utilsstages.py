from rdkit import Chem

__all__ = ("set_atom_pdb_info", )

def set_atom_pdb_info(mol: Chem.Mol, resname: str = "LIG") -> Chem.Mol:
    mol.SetProp("_Name", resname)
    mi = Chem.AtomPDBResidueInfo()
    mi.SetResidueName(resname)
    mi.SetResidueNumber(1)
    mi.SetOccupancy(0.0)
    mi.SetTempFactor(0.0)
    [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
    return mol