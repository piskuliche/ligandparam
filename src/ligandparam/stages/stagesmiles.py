from pathlib import Path
from typing import Union
from typing_extensions import override

from rdkit import Chem

from ligandparam import AbstractStage


class StageSmilesToPDB(AbstractStage):

    @override
    def __init__(self, stage_name: str, input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, input, cwd, *args, **kwargs)
        self.in_smiles = input
        self.out_pdb = Path(kwargs["out_pdb"])
        self.resname = kwargs.get("resname", "LIG")
        self.reduce = kwargs.get("reduce", True)
        self.add_conect = kwargs.get("add_conect", True)
        self.random_seed = kwargs.get("random_seed", None)

    def _execute(self, dry_run=False):
        # First, create the molecule
        mol = Chem.MolFromSmiles(self.smiles)
        if self.reduce:
            mol = Chem.rdmolops.AddHs(mol)
        # All the atoms have their coordinates set to zero. Come up with some values
        params = Chem.AllChem.ETKDGv3()
        if self.random_seed:
            params.randomSeed = self.random_seed
        Chem.AllChem.EmbedMolecule(mol, params)

        # Set metadata and write away
        mol.SetProp("_Name", "LIG")
        mi = Chem.AtomPDBResidueInfo()
        mi.SetResidueName("LIG")
        mi.SetResidueNumber(1)
        mi.SetOccupancy(0.0)
        mi.SetTempFactor(0.0)
        [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
        flavor = 0 if self.add_conect else 2
        #
        Chem.MolToPDBFile(mol, self.out_pdb, flavor=flavor)

        Chem.rdmolfiles.MolToPDBFile(mol, self.out_pdb)

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        raise NotImplementedError

    def _clean(self):
        raise NotImplementedError
