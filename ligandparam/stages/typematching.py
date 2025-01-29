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

    def __init__(self, stage_name: str, name: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, name, cwd, *args, **kwargs)
        
        """ This class updates either (or both) the atom types and names in a mol2 file to match another mol2 file.
        
        Parameters
        ----------
        name : str
            The name of the stage
        in_mol2 : str
            The original mol2 file
        to_update : str
            The file to update
        out_mol2 : str
            The new mol2 file
        update_names : bool
            If True, update the atom names
        update_types : bool
            If True, update the atom types
        update_resname : bool
            If True, update the residue names (only if another update is requested)
        inputoptions : dict
            The input options for the stage
        """
        for opt in ("in_mol2", "out_mol2", "to_update"):
            try:
                setattr(self, opt, kwargs[opt])
            except KeyError:
                raise ValueError(f"ERROR: Please provide {opt} option as a keyword argument.")

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
                        self.logger.debug(f"Atom with {orig_atom.name} has type {orig_atom.type} and will be updated to {new_atom.type}")
                        new_atom.type = orig_atom.type
                if orig_atom.name != new_atom.name:
                    if self.update_names:
                        self.logger.debug(f"Atom with {new_atom.name} will be updated to {orig_atom.name}")
                        new_atom.name = orig_atom.name

            types_mol2 = Path(self.cwd, f"{self.name.stem}.types.mol2")
            if not dry_run:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    Mol2Writer(unew, filename=types_mol2).write()

            ante = Antechamber(cwd=self.cwd, logger=self.logger)
            ante.call(i=types_mol2, fi='mol2',
                      o=self.out_mol2, fo='mol2',
                      pf='y', at=self.atom_type,
                      dry_run = dry_run)


    def _clean(self):
        raise NotImplementedError("clean method not implemented")