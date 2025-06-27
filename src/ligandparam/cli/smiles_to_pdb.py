import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_pdb(smiles, pdb_filename, resname="LIG"):
    # Generate molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        raise RuntimeError("3D coordinate generation failed.")

    # Minimize with MMFF or UFF
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        AllChem.UFFOptimizeMolecule(mol)
    # Set residue name for all atoms
    for atom in mol.GetAtoms():
        atom.SetProp("resName", resname)
    # Write to PDB
    with open(pdb_filename, 'w') as f:
        f.write(Chem.MolToPDBBlock(mol))

def main():
    parser = argparse.ArgumentParser(description="Convert SMILES to PDB with 3D coordinates and minimization.")
    parser.add_argument("-s", "--smiles", required=True, help="Input SMILES string")
    parser.add_argument("-o", "--output", required=True, help="Output PDB filename")
    parser.add_argument("-rn", "--resname", default="LIG", help="Residue name for the ligand (default: LIG)")
    args = parser.parse_args()

    smiles_to_pdb(args.smiles, args.output, args.resname)

if __name__ == "__main__":
    main()