from pathlib import Path
from typing import Optional, Union, Any

from ligandparam.stages import AbstractStage
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.AllChem import ETKDGv3, EmbedMolecule, AlignMol
from typing_extensions import override

from ligandparam.stages import set_atom_pdb_info

class StageSmilesToPDB(AbstractStage):

    @override
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_smiles = main_input
        self.out_pdb = Path(kwargs["out_pdb"])
        self.resname = kwargs.get("resname", "LIG")
        self.reduce = kwargs.get("reduce", True)
        self.add_conect = kwargs.get("add_conect", True)
        self.random_seed = kwargs.get("random_seed", None)

        try:
            self.reference_pdb = Path(kwargs["reference_pdb"]).resolve()
            self.add_required(self.reference_pdb)
            self.normalize_atom_names = True
            self.align = kwargs.get("align", False)
        except KeyError:
            self.normalize_atom_names = False
            self.align = False

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        super()._setup_execution(dry_run=dry_run, nproc=nproc, mem=mem)
        # First, create the molecule
        try:
            mol = Chem.MolFromSmiles(self.in_smiles)
        except Exception as e:
            err_msg = f"Failed to generate an rdkit molecule from input SMILES {self.in_smiles}. Got exception: {e}"
            self.logger.error(err_msg)
            raise RuntimeError(err_msg)

        if self.reduce:
            mol = Chem.rdmolops.AddHs(mol)
        # All the atoms have their coordinates set to zero. Come up with some values
        params = ETKDGv3()
        if self.random_seed:
            params.randomSeed = self.random_seed
        EmbedMolecule(mol, params)

        # Set metadata
        mol = set_atom_pdb_info(mol, self.resname)

        # Normalize the molecule to match the reference PDB
        if self.normalize_atom_names:
            mol = self.normalize_to_reference(mol, self.reference_pdb, self.align)

        flavor = 0 if self.add_conect else 2
        self.logger.info(f"Writing {self.in_smiles} to {self.out_pdb}")

        try:
            Chem.MolToPDBFile(mol, self.out_pdb, flavor=flavor)
        except Exception as e:
            self.logger.error(
                f"Failed to write to  {self.out_pdb}. Got exception: {e}")

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        raise NotImplementedError

    def _clean(self):
        raise NotImplementedError

    def normalize_to_reference(self, mol: Chem.Mol, reference_pdb: Path, align: bool = False) -> Chem.Mol:
        # Normalize the atom names to match the reference PDB
        ref_mol = Chem.MolFromPDBFile(str(reference_pdb))
        if not ref_mol:
            raise ValueError(f"Failed to read reference PDB file {reference_pdb}")

        mcs_mol = self.get_mcs_mol(ref_mol, mol)
        mol_substructure = mol.GetSubstructMatch(mcs_mol)
        ref_substructure = ref_mol.GetSubstructMatch(mcs_mol)
        assert len(mol_substructure) == len(ref_substructure), \
            f"Mismatch in number of common atoms. {len(mol_substructure)} vs {len(ref_substructure)}. This is likely a bug."

        # Get the mapping of the atoms in the target molecule to the reference molecule, and copy names from reference
        atom_map = list(zip(mol_substructure, ref_substructure))
        dict_atom_map = dict(zip(mol_substructure, ref_substructure))
        ref_atoms = list(ref_mol.GetAtoms())
        for a in mol.GetAtoms():
            if a.GetIdx() in dict_atom_map:
                pdb_info = a.GetPDBResidueInfo()
                name = ref_atoms[dict_atom_map[a.GetIdx()]].GetPDBResidueInfo().GetName()
                pdb_info.SetName(name)
                a.SetMonomerInfo(pdb_info)
        if align:
            AlignMol(mol, ref_mol, atomMap=atom_map)
        return mol

    @staticmethod
    def get_mcs_mol(ref_mol: Chem.Mol, mol: Chem.Mol) -> Chem.Mol:
        """ Copy the atom names from the reference molecule to the target molecule. """
        mcs = rdFMCS.FindMCS([ref_mol, mol])
        common_mol = Chem.rdmolfiles.MolFromSmarts(mcs.smartsString)
        return common_mol
