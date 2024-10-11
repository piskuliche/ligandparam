import MDAnalysis as mda
import numpy as np

from ligand_param.stages.abstractstage import AbstractStage

class StageNormalizeCharges(AbstractStage):

    def __init__(self, name, mol2file=None, netcharge=None) -> None:
        self.name = name
        if mol2file is None:
            raise ValueError("Error (Stage {self.name}): mol2 file not set")
        self.mol2file = mol2file

        if netcharge is None:
            print("Net charge not set. Defaulting to 0.0")
            self.netcharge = 0.0
        else:
            self.netcharge = netcharge

        return
    

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage


    def _execute(self, dry_run=False):
        print(f"Executing {self.name} with netcharge={self.netcharge}")
        u = mda.Universe(self.mol2file)
        total_charge = sum(u.atoms.charges)
        charge_diff = total_charge - self.netcharge
        print("Total charge: ", total_charge)
        if np.abs(charge_diff) > 0.0:
            print("Normalizing Charges")
            newcharge = self.normalize_charges(u.atoms.charges, self.netcharge)
            u.atoms.charges = newcharge

        total_charge = sum(u.atoms.charges)
        print("Total charge after normalization: ", total_charge)

    def _clean(self):
        raise NotImplementedError("clean method not implemented")


    def normalize_charges(self, charges, netcharge):
        raise NotImplementedError("normalize_charges method not implemented")
        