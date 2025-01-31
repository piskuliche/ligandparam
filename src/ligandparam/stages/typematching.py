import warnings
from typing import Union

import numpy as np
import MDAnalysis as mda

from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.io.coordinates import Mol2Writer
from ligandparam.log import get_logger


class StageUpdate(AbstractStage):
    """ This class updates either (or both) the atom types and names in a mol2 file to match another mol2 file.

    Parameters
    ----------

    """

    def __init__(self, stage_name: str, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, in_filename, cwd, *args, **kwargs)
        self.in_mol2 = Path(in_filename)
        self.out_mol2 = Path(kwargs["out_mol2"])

        self.to_update = Path(kwargs["to_update"])
        self.update_names = kwargs.get('update_names', False)
        self.update_types = kwargs.get('update_types', False)
        self.update_resname = kwargs.get('update_resname', False)
        self.atom_type = kwargs.get('atom_type', 'gaff2')

        self.add_required(Path(self.in_mol2))
        self.add_required(Path(self.to_update))

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def _execute(self, dry_run=False) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            if not self.update_names and not self.update_types:
                self.logger.debug("No updates requested. Exiting.")
                return
            if self.update_names and self.update_types:
                self.logger.debug("Both updates requested. This will update both atom names and types.")
            elif self.update_names:
                self.logger.debug("Only updating atom names.")
            elif self.update_types:
                self.logger.debug("Only updating atom types.")

            uorig = mda.Universe(self.in_mol2, format="mol2")
            unew = mda.Universe(self.to_update, format="mol2")
            if self.update_resname:
                unew.residues.resnames = uorig.residues.resnames
            for orig_atom, new_atom in zip(uorig.atoms, unew.atoms):
                if orig_atom.type != new_atom.type:
                    if self.update_types:
                        self.logger.debug(
                            f"Atom with {orig_atom.name} has type {orig_atom.type} and will be updated to {new_atom.type}")
                        new_atom.type = orig_atom.type
                if orig_atom.name != new_atom.name:
                    if self.update_names:
                        self.logger.debug(f"Atom with {new_atom.name} will be updated to {orig_atom.name}")
                        new_atom.name = orig_atom.name

            types_mol2 = Path(self.cwd, f"{self.out_mol2.stem}.types.mol2")
            if not dry_run:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    Mol2Writer(unew, filename=types_mol2).write()

            ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
            ante.call(i=types_mol2, fi='mol2',
                      o=self.out_mol2, fo='mol2',
                      pf='y', at=self.atom_type,
                      gn=f"%nproc={self.nproc}", gm=f"%mem={self.mem}MB",
                      dry_run=dry_run)

    def _clean(self):
        raise NotImplementedError("clean method not implemented")
