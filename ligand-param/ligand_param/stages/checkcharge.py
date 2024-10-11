import numpy as np
import MDAnalysis as mda

from ligand_param.stages.abstractstage import AbstractStage

class StageCheckCharge(AbstractStage):
    """ This is an abstract class for all the stages. """
    def __init__(self, name, filename = None, filetype = None, netcharge=0.0) -> None:
        self.name = name
        self.filename = filename
        self.filetype = filetype
        self.netcharge = netcharge
        pass
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def _execute(self, dry_run=False):
        u = mda.Universe(self.filename, format=self.filetype)
        total_charge = sum(u.atoms.charges)
        print("Total charge: ", total_charge)
        compare_charge = self.netcharge - total_charge
        if np.round(compare_charge, 2) != 0.0:
            print("Error: Total charge does not match the net charge")
            print(f"Net charge: {self.netcharge}, Total charge: {total_charge}")
            print(f"Charge difference: {compare_charge}")
            raise ValueError("Error: Total charge does not match the net charge")
        return

    def _clean(self):
        raise NotImplementedError("clean method not implemented")

