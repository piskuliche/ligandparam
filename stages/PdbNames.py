from pathlib import Path
from typing import Optional, Union, Any

from ligandparam.stages.AbstractStage import AbstractStage
from rdkit import Chem
from typing_extensions import override

from ligandparam.io.Smiles import normalize_to_reference as normalize_mol_to_reference
from ligandparam.stages.StageUtils import set_atom_pdb_info


class StagePdbNameFixer(AbstractStage):
    """Fix PDB atom names and related metadata in ligand files."""

    @override
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_pdb = main_input
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

    def _run(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        try:
            mol = Chem.MolFromPDBFile(str(self.in_pdb), removeHs=False)
        except Exception as e:
            err_msg = (
                f"Failed to generate an rdkit molecule from input PDB "
                f"{self.in_pdb}. Got exception: {e}"
            )
            self.logger.error(err_msg)
            raise RuntimeError(err_msg) from e

        mol = set_atom_pdb_info(mol, self.resname)

        if self.normalize_atom_names:
            mol = normalize_mol_to_reference(
                mol,
                self.reference_pdb,
                align=self.align,
                remove_hs=False,
                logger=self.logger,
            )

        flavor = 0 if self.add_conect else 2
        self.logger.info(f"Writing {self.in_pdb} to {self.out_pdb}")

        try:
            Chem.MolToPDBFile(mol, str(self.out_pdb), flavor=flavor)
        except Exception as e:
            self.logger.error(f"Failed to write to  {self.out_pdb}. Got exception: {e}")


# Back-compat alias (CLI / older recipes).
PDB_Name_Fixer = StagePdbNameFixer
