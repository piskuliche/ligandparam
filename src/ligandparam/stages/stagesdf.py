from pathlib import Path
from typing import Union
from typing_extensions import override

from rdkit import Chem

from ligandparam.stages import AbstractStage


class StageSDFToPDB(AbstractStage):

    @override
    def __init__(self, stage_name: str, input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, input, cwd, *args, **kwargs)
        self.in_sdf = Path(input)
        self.out_pdb = Path(kwargs["out_pdb"])
        self.resname = kwargs.get("resname", "LIG")
        self.removeHs = kwargs.get("removeHs", False)
        self.add_conect = kwargs.get("add_conect", True)

    def _execute(self, dry_run=False):
        # First, create the molecule
        try:
            mols = Chem.SDMolSupplier(self.in_sdf, removeHs=False)
        except Exception as e:
            err_msg = f"Failed to generate an rdkit molecule from input SDF {self.in_sdf}. Got exception: {e}"
            self.logger.error(err_msg)
            raise RuntimeError(err_msg)

        if len(mols) > 1:
            err_msg = f"Input SDF,  {self.in_sdf}, has more than 1 molecule"
            self.logger.error(err_msg)
            raise ValueError(err_msg)
        mol = mols[0]

        # Set metadata and write away
        mol.SetProp("_Name", self.resname)
        mi = Chem.AtomPDBResidueInfo()
        mi.SetResidueName(self.resname)
        mi.SetResidueNumber(1)
        mi.SetOccupancy(0.0)
        mi.SetTempFactor(0.0)
        [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
        flavor = 0 if self.add_conect else 2
        self.logger.info(f"Writing {self.in_sdf} to {self.out_pdb}")

        try:
            Chem.MolToPDBFile(mol, self.out_pdb, flavor=flavor)
        except Exception as e:
            self.logger.error(
                f"Failed to write to  {self.out_pdb}. Got exception: {e}")

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        raise NotImplementedError

    def _clean(self):
        raise NotImplementedError
