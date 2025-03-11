from pathlib import Path
from typing import Union, Iterable

import numpy as np
from typing_extensions import override
from collections import Counter

from rdkit import Chem

from ligandparam.stages import AbstractStage


class SDFToPDBBatch(AbstractStage):

    @override
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_sdf = Path(main_input)

        self.removeHs = kwargs.get("removeHs", False)
        self.add_conect = kwargs.get("add_conect", True)

        self.out_pdb_template = kwargs.get("out_pdb_template", None)
        self.out_pdbs = kwargs.get("out_pdbs", None)
        self.resnames = kwargs.get("resnames", None)
        self.resname = kwargs.get("resname", None)

    def execute(self, dry_run=False, nproc=1, mem=512):
        # First, create the molecule
        try:
            mols = Chem.SDMolSupplier(self.in_sdf, removeHs=False)
        except Exception as e:
            err_msg = f"Failed to generate an rdkit molecule from input SDF {self.in_sdf} Got exception: {e}"
            self.logger.error(err_msg)
            raise RuntimeError(err_msg)

        # Set up names and paths
        if self.resnames is None:
            if self.resname is None:
                self.resnames = [mol.GetProp("_Name")[:3] for mol in mols]
            elif self.resname:
                self.resnames = [self.resname for _ in mols]
        if self.out_pdbs is None:
            if self.out_pdb_template is None:
                filenames = [f'{mol.GetProp("_Name")}.pdb' for mol in mols]
                counts = Counter(filenames)
                if np.all(np.array(list(counts.values())) == 1):
                    self.out_pdbs = [self.cwd / fn for fn in filenames]
                else:
                    err_msg = f"Multiple molecules with the same name in {self.in_sdf} Please provide `out_pdbs` or `out_pdb_template`."
                    self.logger.error(err_msg)
                    raise ValueError(err_msg)
            else:
                out_dir = self.out_pdb_template.parent
                label = self.out_pdb_template.stem
                self.out_pdbs = [out_dir / f"{label}_{i}.pdb" for i in range(0, len(self.resnames))]

        if len(self.resnames) != len(self.out_pdbs) or len(self.resnames) != len(mols):
            err_msg = f"Lengths of `out_pdbs`, `resnames`, and mols don't match: {len(self.out_pdbs)}, {len(self.resnames)}, and {len(mols)}"
            self.logger.error(err_msg)
            raise ValueError(err_msg)

        # Write each mol to a different PDB
        flavor = 0 if self.add_conect else 2
        for mol, pdb, resname in zip(mols, self.out_pdbs, self.resnames):
            # Set metadata and write away
            mol.SetProp("_Name", resname)
            mi = Chem.AtomPDBResidueInfo()
            mi.SetResidueName(resname)
            mi.SetResidueNumber(1)
            mi.SetOccupancy(0.0)
            mi.SetTempFactor(0.0)
            [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
            self.logger.info(f"Writing {self.in_sdf} to {pdb}")

            try:
                Chem.MolToPDBFile(mol, pdb, flavor=flavor)
            except Exception as e:
                self.logger.error(
                    f"Failed to write to  {pdb}. Got exception: {e}")

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        raise NotImplementedError

    def _clean(self):
        raise NotImplementedError
