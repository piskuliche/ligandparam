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

        self.source_mol2 = Path(kwargs["source_mol2"])
        self.update_names = kwargs.get('update_names', False)
        self.update_types = kwargs.get('update_types', False)
        self.update_resname = kwargs.get('update_resname', False)
        self.update_charges = kwargs.get('update_charges', False)
        self.atom_type = kwargs.get('atom_type', 'gaff2')
        if "molname" in kwargs:
            self.additional_args = {"rn": kwargs["molname"]}
        else:
            self.additional_args = {}

        self.add_required(Path(self.in_mol2))
        self.add_required(Path(self.source_mol2))

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def _execute(self, dry_run=False) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            if not (self.update_names or self.update_types or self.update_charges):
                self.logger.debug("No updates requested. Exiting.")
                return
            if self.update_names and self.update_types:
                self.logger.debug("Both updates requested. This will update both atom names and types.")
            elif self.update_names:
                self.logger.debug("Only updating atom names.")
            elif self.update_types:
                self.logger.debug("Only updating atom types.")

            dest_u = mda.Universe(self.in_mol2, format="mol2")
            source_u = mda.Universe(self.source_mol2, format="mol2")
            if self.update_resname:
                dest_u.residues.resnames = source_u.residues.resnames
            for orig_atom, new_atom in zip(source_u.atoms, dest_u.atoms):
                if self.update_types:
                    self.logger.debug(
                        f"Atom {orig_atom.name} with type {orig_atom.type} and will be updated to {new_atom.type}")
                    new_atom.type = orig_atom.type
                if self.update_names:
                    self.logger.debug(f"Atom {new_atom.name} will named {orig_atom.name}")
                    new_atom.name = orig_atom.name
                if self.update_charges:
                    self.logger.debug(f"Atom {new_atom.name}'s charge ({orig_atom.charge}) will be updated to {orig_atom.name}")
                    new_atom.charge = orig_atom.charge

            types_mol2 = Path(self.cwd, f"{self.out_mol2.stem}.types.mol2")
            if not dry_run:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    Mol2Writer(dest_u, filename=types_mol2).write()

            ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
            ante.call(i=types_mol2, fi='mol2',
                      o=self.out_mol2, fo='mol2',
                      pf='y', at=self.atom_type,
                      gn=f"%nproc={self.nproc}", gm=f"%mem={self.mem}MB",
                      dry_run=dry_run, **self.additional_args)

    def _clean(self):
        raise NotImplementedError("clean method not implemented")
