from pathlib import Path
from typing import Optional,  Union, Any
from typing_extensions import override

from rdkit import Chem

from ligandparam.stages import AbstractStage


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

    def execute(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> Any:
        # First, create the molecule
        try:
            mol = Chem.MolFromSmiles(self.in_smiles)
        except Exception as e:
            err_msg = f"Failed to generate an rdkit molecule from input SDF {self.in_sdf}. Got exception: {e}"
            self.logger.error(err_msg)
            raise RuntimeError(err_msg)

        if self.reduce:
            mol = Chem.rdmolops.AddHs(mol)
        # All the atoms have their coordinates set to zero. Come up with some values
        params = Chem.AllChem.ETKDGv3()
        if self.random_seed:
            params.randomSeed = self.random_seed
        Chem.AllChem.EmbedMolecule(mol, params)

        # Set metadata and write away
        mol.SetProp("_Name", self.resname)
        mi = Chem.AtomPDBResidueInfo()
        mi.SetResidueName(self.resname)
        mi.SetResidueNumber(1)
        mi.SetOccupancy(0.0)
        mi.SetTempFactor(0.0)
        [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
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
