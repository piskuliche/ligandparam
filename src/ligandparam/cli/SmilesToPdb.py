"""CLI to convert a SMILES string to a 3D PDB."""

from __future__ import annotations

import argparse

from ligandparam.io.Smiles import PDBFromSMILES


def smiles_to_pdb(smiles: str, pdb_filename: str, resname: str = "LIG") -> None:
    """Embed a SMILES string and write a cleaned PDB (shared ``io.smiles`` path)."""
    builder = PDBFromSMILES(resname, smiles)
    builder.mol_from_smiles(addHs=True)
    builder.write_pdb(pdb_filename, minimize=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert SMILES to PDB with 3D coordinates."
    )
    parser.add_argument("-s", "--smiles", required=True, help="Input SMILES string")
    parser.add_argument("-o", "--output", required=True, help="Output PDB filename")
    parser.add_argument(
        "-rn",
        "--resname",
        default="LIG",
        help="Residue name for the ligand (default: LIG)",
    )
    args = parser.parse_args(argv)
    smiles_to_pdb(args.smiles, args.output, args.resname)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
