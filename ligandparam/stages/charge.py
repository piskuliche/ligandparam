import shutil

import numpy as np
import MDAnalysis as mda

from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.io.coordinates import Mol2Writer

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

        self.orig_mol2 = orig_mol2
        self.new_mol2 = new_mol2

        if self.orig_mol2 == self.new_mol2:
            raise ValueError("ERROR: Original and new mol2 files are the same. Please provide different files.")
        
        if charge_source is "multistate":
            self.charge_source = "respfit.out"
            self.charge_column = 3
        else:
            if charge_source is not None:
                self.charge_source = charge_source
            else:
                raise ValueError("ERROR: Please provide a charge source file.")
            
            if charge_column is not None:
                self.charge_column = charge_column
            else:
                raise ValueError("ERROR: Please provide a charge column.")
        
        self.add_required(self.orig_mol2)
        self.add_required(self.charge_source)

        return
    
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
                raise ValueError("Error: Number of charges does not match the number of atoms.")
            u.atoms.charges = charges
            # Write the Mol2 temporary file
            Mol2Writer(u, self.base_cls.base_name + ".tmpresp.mol2", selection="all").write()
        
        ante = Antechamber()
        ante.call(i=self.base_cls.base_name + ".tmpresp.mol2", fi='mol2',
                  o=self.new_mol2, fo='mol2',
                  pf='y', at=self.base_cls.atom_type,
                  dry_run = dry_run)

        return

    def _clean(self):
        raise NotImplementedError("clean method not implemented")

    
class StageNormalizeCharge(AbstractStage):
    """ This class normalizes the charges to the net charge. 
    
    This class works by calculating the charge difference, and then normalizing the charges
    based on the overall precision that you select, by adjusting each atom charge by the precision
    until the charge difference is zero."""

    def __init__(self, name, base_cls=None, orig_mol2=None, new_mol2=None, precision=0.0001):
        """ Initialize the StageNormalizeCharge class.
        
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
        precision : float
            The precision of the charge normalization

        Raises
        ------
        ValueError
            If the original mol2 file is not provided
        ValueError
            If the new mol2 file is not provided

        """
        self.name = name
        self.base_cls = base_cls
        if orig_mol2 is not None:
            self.orig_mol2 = orig_mol2
        else:
            raise ValueError("Please provide an original filename.")
        
        if new_mol2 is not None:
            self.new_mol2 = new_mol2
        else:
            raise ValueError("Please provide a new filename.")
        

        self.precision = precision
        self.decimals = len(str(precision).split(".")[1])

        self.add_required(orig_mol2)


    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def _execute(self, dry_run=False):
        """ Execute the stage. 

        Raises
        ------
        ValueError
            If the charge normalization fails
        
        TODO: Check what happens when netcharge is nonzero
        TODO: Check what happens when charge difference is larger than the number of atoms

        """
        import warnings
        # Supress the inevitable mol2 file warnings.
        warnings.filterwarnings("ignore")
        print("-> Checking charges")
        print(f"-> Normalizing charges to {self.base_cls.net_charge}")
        print(f"-> Precision {self.precision} with {self.decimals} decimals")

        u = mda.Universe(self.orig_mol2, format='mol2')
        total_charge, charge_difference = self.check_charge(u.atoms.charges)
        
        if charge_difference != 0.0:
            print("-> Normalizing charges")
            charges = self.normalize(u.atoms.charges, charge_difference)
            new_total, new_diff = self.check_charge(charges)
            if new_diff != 0.0:
                raise ValueError("Error: Charge normalization failed.")
            else:
                u.atoms.charges = charges
        else:
            print("-> Charges are already normalized")
            return
        if not dry_run:
            Mol2Writer(u, self.base_cls.base_name + ".tmpnorm.mol2", selection="all").write()

        ante = Antechamber()
        ante.call(i=self.base_cls.base_name + ".tmpnorm.mol2", fi='mol2',
                  o=self.new_mol2, fo='mol2',
                  pf='y', at=self.base_cls.atom_type,
                  dry_run = dry_run)

    def _clean(self):
        raise NotImplementedError("clean method not implemented")

    def normalize(self, charges, charge_difference):
        """ This function normalizes the charges to the net charge. 
        
        Parameters
        ----------
        charges : np.array
            The charges
        charge_difference : float
            The charge difference
        
        Returns
        -------
        charges : np.array
            The normalized charges
        """

        count = np.round(np.abs(charge_difference)/self.precision)
        adjust = np.round(charge_difference/count, self.decimals)
        natoms = len(charges)
        # Choosing charges closest to zero.
        sorted_indices = np.argsort(np.abs(charges))
        # Flip the order to choose the largest charges first.
        sorted_indices = sorted_indices[::-1]
        for i in range(int(count)):
            atom_idx = i % natoms
            charges[sorted_indices[atom_idx]] += adjust
        return charges

    def check_charge(self, charges):
        """ This function checks the total charge and the charge difference.
        
        Parameters
        ----------
        charges : np.array
            The charges
        
        Returns
        -------
        total_charge : float
            The total charge
        charge_difference : float
            The charge difference
        """
        total_charge = np.round(sum(charges), self.decimals)
        charge_difference = np.round(self.base_cls.net_charge - total_charge, self.decimals)
        print(f"-> Total Charge is {total_charge}")
        print(f"-> Charge difference is {charge_difference}")
        return total_charge, charge_difference




