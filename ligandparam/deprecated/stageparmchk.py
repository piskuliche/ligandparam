import MDAnalysis as mda
import numpy as np

from ligandparam.abstractstage import AbstractStage
from ligandparam.interfaces import ParmChk

class StageParmChk(AbstractStage):

    def __init__(self, name, base_cls=None) -> None:
        self.name = name
        self.base_cls = base_cls
        return
    

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage


    def execute(self, dry_run=False):
        print(f"Executing {self.name} with netcharge={self.base_cls.net_charge}")
        parm = ParmChk()
        parm.call(i=self.base_cls.name+'.resp.mol2', f="mol2",
                  o=self.base_cls.name+'.frcmod', 
                  s=2, dry_run = dry_run)
        return
