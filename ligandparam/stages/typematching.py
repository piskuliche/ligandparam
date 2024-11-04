import numpy as np
import MDAnalysis as mda

from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.io.coordinates import Mol2Writer

class StageUpdate(AbstractStage):
    
    def __init__(self, name, orig_mol2=None, to_update=None, new_mol2=None, update_names=False, update_types=False, update_resname=False, inputoptions=None) -> None:
        """ This class updates either (or both) the atom types and names in a mol2 file to match another mol2 file.
        
        Parameters
        ----------
        name : str
            The name of the stage
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
        update_resname : bool
            If True, update the residue names (only if another update is requested)
        inputoptions : dict
            The input options for the stage
        """
        self.name = name
        self._parse_inputoptions(inputoptions)
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
        self.update_resname = update_resname

        self._add_required(orig_mol2)
        self._add_required(to_update)
        self._add_output(new_mol2)

    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def _execute(self, dry_run=False):
        """ Execute the update stage, which updates the atom names and types in a mol2 file.
        
        This stage will update atom names/atom types/resnames in a mol2 file to match another mol2 file, if requested. These can be updated individually or together.
        
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would"""
        import warnings
        # Supress the inevitable mol2 file warnings.
        warnings.filterwarnings("ignore")

        if not self.update_names and not self.update_types:
            print("No updates requested. Exiting.")
            return
        if self.update_names and self.update_types:
            print("Both updates requested. This will update both atom names and types.")
        elif self.update_names:
            print("Only updating atom names.")
        elif self.update_types:
            print("Only updating atom types.")
        

        uorig = mda.Universe(self.orig_mol2, format="mol2")
        unew = mda.Universe(self.to_update, format="mol2")
        if self.update_resname:
            unew.residues.resnames = uorig.residues.resnames
        for orig_atom, new_atom in zip(uorig.atoms, unew.atoms):
            if orig_atom.type != new_atom.type:
                if self.update_types:
                    print(f"Atom with {orig_atom.name} has type {orig_atom.type} and will be updated to {new_atom.type}")
                    new_atom.type = orig_atom.type
            if orig_atom.name != new_atom.name:
                if self.update_names:
                    print(f"Atom with {new_atom.name} will be updated to {orig_atom.name}")
                    new_atom.name = orig_atom.name
            
        if not dry_run:
            Mol2Writer(unew, filename=f"{self.base_name}.types.mol2").write()

        ante = Antechamber()
        ante.call(i=self.base_name + ".types.mol2", fi='mol2',
                  o=self.new_mol2, fo='mol2',
                  pf='y', at=self.atom_type,
                  dry_run = dry_run)
        return

    def _clean(self):
        raise NotImplementedError("clean method not implemented")