
from abc import abstractmethod
from abstractstage import AbstractStage
from antechamber_interface import Antechamber

class StageInitialize(AbstractStage):
    def __init__(self, name, base_name=None, net_charge=None, atom_type=None) -> None:
        self.name = name

        if base_name is None:
            raise ValueError("Error (Stage {self.name}): Base name not set")
        
        self.base_name = base_name
        if net_charge is None:
            print("Net charge not set. Defaulting to 0.0")
            self.net_charge = 0.0
        else:
            self.net_charge = net_charge
        
        if atom_type is None:
            print("Atom type not set. Defaulting to gaff2")
            self.atom_type = 'gaff2'
        else:
            self.atom_type = atom_type
        
        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False):
        ante = Antechamber()
        ante.call(i=self.base_name+'.pdb', fi='pdb',
                  o=self.base_name+'.mol2', fo='mol2',
                  c='bcc', nc=self.net_charge,
                  pf='y', at=self.atom_type,
                  run=(not dry_run))