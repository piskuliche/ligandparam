import shutil

import numpy as np
import MDAnalysis as mda

from pathlib import Path

from ligand_param.stages.abstractstage import AbstractStage
from ligand_param.interfaces import Antechamber

class StageUpdateCharge(AbstractStage):
    """ This class creates a new mol2 file with updated charges. """
    def __init__(self, name, base_cls=None, orig_mol2=None, new_mol2=None, charge_source="multistate", charge_column=None) -> None:
        """ Initialize the StageUpdateCharge class.
        
        Parameters
        ----------
        name : str
            The name of the stage
        base_cls : Ligand
            The base class of the ligand
        orig_mol2 : str
            The original mol2 file
        new_mol2 : str
            The new mol2 file
        charge_source : str
            The source of the charges
        charge_column : int
            The column of the charges
        """
        self.name = name
        self.base_cls = base_cls

        if orig_mol2 is not None:
            self.orig_mol2 = orig_mol2
        else:
            self.orig_mol2 = self.base_cls.base_name + ".mol2"
        
        if new_mol2 is not None:
            self.new_mol2 = new_mol2
        else:
            self.new_mol2 = self.base_cls.base_name + ".resp.mol2"

        if self.orig_mol2 == self.new_mol2:
            raise ValueError("Original and new mol2 files are the same. Please provide different files.")
        if charge_source is "multistate":
            self.charge_source = "respfit.out"
            self.charge_column = 3
        else:
            if charge_source is not None:
                self.charge_source = charge_source
            else:
                raise ValueError("Please provide a charge source file.")
            
            if charge_colunm is not None:
                self.charge_column = charge_column
            else:
                raise ValueError("Please provide a charge column.")
        pass
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def _execute(self, dry_run=False):
        import warnings
        # Supress the inevitable mol2 file warnings.
        warnings.filterwarnings("ignore")
        if Path(self.charge_source).exists():
            charges = np.genfromtxt(self.charge_source, usecols=(self.charge_column), unpack=True)
        else:
            raise FileNotFoundError(f"File {self.charge_source} not found.")

        if not dry_run:
            u = mda.Universe(self.orig_mol2, format='mol2')
            if len(charges) != len(u.atoms):
                raise ValueError("Number of charges does not match the number of atoms.")
            u.atoms.charges = charges
            ag = u.select_atoms("all")
            ag.write(self.base_cls.base_name + ".tmpresp.mol2")
        
        ante = Antechamber()
        ante.call(i=self.base_cls.base_name + ".tmpresp.mol2", fi='mol2',
                  o=self.new_mol2, fo='mol2',
                  pf='y', at=self.base_cls.atom_type,
                  dry_run = dry_run)

        return

    def _clean(self):
        raise NotImplementedError("clean method not implemented")