import numpy as np
import MDAnalysis as mda

from pathlib import Path

from ligand_param.stages.abstractstage import AbstractStage
from ligand_param.interfaces import Antechamber
from ligand_param.io.coordinates import Mol2Writer

class StageUpdate(AbstractStage):
    """" This class updates either (or both) the atom types and names in a mol2 file to match another mol2 file. """
    def __init__(self, name, base_cls=None, orig_mol2=None, to_update=None, new_mol2=None, update_names=False, update_types=False) -> None:
        """ Initialize the StageUpdate class.
        
        Parameters
        ----------
        name : str
            The name of the stage
        base_cls : Ligand
            The base class of the ligand
        orig_mol2 : str
            The original mol2 file
        to_update : str
            The file to update
        new_mol2 : str
            The new mol2 file
        update_names : bool
            If True, update the atom names
        update_types : bool
            If True, update the atom types
        """
        self.name = name
        self.base_cls = base_cls
        if orig_mol2 is not None:
            self.orig_mol2 = orig_mol2
        else:
            raise ValueError("Please provide the original types file.")

        if to_update is not None:
            self.to_update = to_update
        else:
            raise ValueError("Please provide the file to update.")
        
        if new_mol2 is not None:
            self.new_mol2 = new_mol2
        else:
            raise ValueError("Please provide the new types file.")
        
        self.update_names = update_names
        self.update_types = update_types

    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def _execute(self, dry_run=False):
        import warnings
        # Supress the inevitable mol2 file warnings.
        warnings.filterwarnings("ignore")

        if not self.update_names and not self.update_types:
            print("No updates requested. Exiting.")
            return

        uorig = mda.Universe(self.orig_mol2, format="mol2")
        unew = mda.Universe(self.to_update, format="mol2")
        for orig_atom, new_atom in zip(uorig.atoms, unew.atoms):
            if orig_atom.type != new_atom.type:
                if self.update_types:
                    print(f"Atom {orig_atom.name} has type {orig_atom.type} and will be updated to {new_atom.type}")
                    new_atom.type = orig_atom.type
            if orig_atom.type != new_atom.type:
                if self.update_names:
                    print(f"Atom with {new_atom.name} will be updated to {orig_atom.type}")
                    new_atom.name = orig_atom.name
            
        if not dry_run:
            Mol2Writer(unew, filename=f"{self.base_cls.base_name}.types.mol2").write()

        ante = Antechamber()
        ante.call(i=self.base_cls.base_name + ".types.mol2", fi='mol2',
                  o=self.new_mol2, fo='mol2',
                  pf='y', at=self.base_cls.atom_type,
                  dry_run = dry_run)
        return

    def _clean(self):
        raise NotImplementedError("clean method not implemented")