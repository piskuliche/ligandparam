from typing import Optional, List

from ligandparam.stages import AbstractStage
from ligandparam.io.leapIO import LeapWriter
from ligandparam.interfaces import Leap
from ligandparam.log import get_logger


class StageLeap(AbstractStage):
    
        def __init__(self, stage_name, base_cls=None, leaprc: Optional[List] = None) -> None:
            self.stage_name = stage_name
            self.base_cls = base_cls
            if leaprc is None:
                self.leaprc = []
            else:
                self.leaprc = leaprc

            return
        
    
        def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
            return stage
        
    
        def execute(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> Any:
            print(f"Executing {self.stage_name} with netcharge={self.base_cls.net_charge}")
            leapgen = LeapWriter("param")
            for rc in self.leaprc:
                leapgen.add_leaprc(rc)
            leapgen.add_line(f"loadamberparams {self.base_cls.name}.frcmod")
            leapgen.add_line(f"mol = loadmol2 {self.base_cls.name}.resp.mol2")
            leapgen.add_line(f"saveOff mol {self.base_cls.name}.lib")
            leapgen.add_line("quit")
            leapgen.write()

            leap = Leap()
            leap.call(f="tleap.param.in", dry_run = dry_run)
            return
